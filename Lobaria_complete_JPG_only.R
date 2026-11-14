#!/usr/bin/env Rscript
# Complete Lobaria analysis and NUMBERED manuscript figures: JPG ONLY.
# Figure1–Figure5 and Supp-F1–Supp-F21; all plot modules embedded.
# Place in Lobaria_reanalysis and set that project as your working directory.
# source("Lobaria_complete_JPG_only.R")
# Requires the input folder. python_results is optional for common draws/checks.
# RUN_MODE="auto": use completed native R tables/objects when available; otherwise full analysis.
# RUN_MODE="full": rerun all statistics. RUN_MODE="figures": require completed results.
# No external R scripts are sourced. Existing PDF/PNG files are neither rewritten nor deleted.
# Exports: R_results/Final_Manuscript_JPG/Main/Figure1.jpg ... Figure5.jpg
#          R_results/Final_Manuscript_JPG/Supplementary/Supp-F1.jpg ... Supp-F21.jpg
# Captions and numbering index accompany 600-dpi JPG exports.
# Final physical sizes mostly 5.2–7.2 inches wide; inspect the exported JPGs before submission.
# Palette and inferential cautions are documented in figure captions.
# ASM figure guidance: https://journals.asm.org/writing-your-paper?journal=spectrum
RUN_MODE <- getOption("lobaria.run_mode", "auto")
# Publication figure suite. Embedded in the complete workflow.
publication_figures <- function(context=parent.frame()) evalq({
 suppressPackageStartupMessages({library(ggplot2);library(ggrepel)})
 fdir <- file.path(out,'Final_Manuscript_JPG');dir.create(fdir,showWarnings=FALSE);dir.create(file.path(fdir,'Main'),showWarnings=FALSE);dir.create(file.path(fdir,'Supplementary'),showWarnings=FALSE)
 popcols <- setNames(c('#0072B2','#56B4E9','#009E73','#E69F00','#D55E00','#CC79A7','#332288'),pops)
 popshape <- setNames(c(16,17,15,18,8,3,4),pops)
 partner_names <- c(Cyanobiont='Cyanobiont',Mycobiont='Mycobiont',Noncyano_bacteria='Non-cyanobacterial bacteria',Phycobiont='Phycobiont')
 partnercols <- c(Cyanobiont='#0072B2',Mycobiont='#882255',Noncyano_bacteria='#D55E00',Phycobiont='#009E73')
 theme_set(theme_classic(base_size=11,base_family='sans')+theme(axis.title=element_text(size=12),axis.text=element_text(size=10,color='black'),plot.title=element_text(size=12,face='bold'),plot.subtitle=element_text(size=9),legend.title=element_text(size=10),legend.text=element_text(size=9),strip.text=element_text(size=11,face='bold'),strip.background=element_blank(),legend.position='bottom',plot.margin=margin(8,10,8,8)))
 # Arrange ggplot grobs with base grid: no patchwork operators or version coupling.
 combine_plots <- function(plots,nrow=NULL,ncol=NULL,guides=NULL,legend_position='bottom') {
  if(inherits(plots,'ggplot'))plots<-list(plots)
  if(is.null(ncol))ncol<-if(is.null(nrow))1L else ceiling(length(plots)/nrow)
  if(is.null(nrow))nrow<-ceiling(length(plots)/ncol)
  legend<-NULL
  if(identical(guides,'collect')){
   first<-ggplot2::ggplotGrob(plots[[1]]);ii<-grep('^guide-box',first$layout$name)
   ii<-ii[vapply(first$grobs[ii],function(x)inherits(x,'gtable'),logical(1))]
   if(length(ii))legend<-first$grobs[[ii[1]]]
   plots<-lapply(plots,function(p)p+ggplot2::theme(legend.position='none'))
  }
  grobs<-lapply(plots,function(p)if(inherits(p,'grob'))p else ggplot2::ggplotGrob(p))
  if(all(vapply(grobs,function(x)!is.null(x$widths),logical(1))) && length(unique(vapply(grobs,function(x)length(x$widths),integer(1))))==1L){
   common<-do.call(grid::unit.pmax,lapply(grobs,function(x)x$widths))
   grobs<-lapply(grobs,function(x){x$widths<-common;x})
  }
  children<-lapply(seq_along(grobs),function(i)grid::grobTree(grobs[[i]],vp=grid::viewport(layout.pos.row=(i-1L)%/%ncol+1L,layout.pos.col=(i-1L)%%ncol+1L)))
  body<-grid::gTree(children=do.call(grid::gList,children),vp=grid::viewport(layout=grid::grid.layout(nrow,ncol)))
  if(is.null(legend))return(body)
  if(legend_position=='right'){
   return(grid::gTree(children=grid::gList(grid::grobTree(body,vp=grid::viewport(layout.pos.col=1)),grid::grobTree(legend,vp=grid::viewport(layout.pos.col=2))),vp=grid::viewport(layout=grid::grid.layout(1,2,widths=grid::unit.c(grid::unit(1,'null'),grid::grobWidth(legend)+grid::unit(.12,'in'))))))
  }
  grid::gTree(children=grid::gList(grid::grobTree(body,vp=grid::viewport(layout.pos.row=1)),grid::grobTree(legend,vp=grid::viewport(layout.pos.row=2))),vp=grid::viewport(layout=grid::grid.layout(2,1,heights=grid::unit.c(grid::unit(1,'null'),grid::grobHeight(legend)+grid::unit(.15,'in')))))
 }

 # Build plots first; only the final numbered figures are exported below.
 plotcache <- list()
 savefig <- function(p,name,caption,width=7.2,height=5.4) {
  plotcache[[name]] <<- list(plot=p,caption=caption,width=width,height=height)
 }
 wrap <- function(x,n=35)vapply(x,function(y)paste(strwrap(y,n),collapse='\n'),character(1))
 long <- function(x) {d<-as.data.frame(as.table(x),stringsAsFactors=FALSE);names(d)<-c('sample','feature','value');d$population<-factor(as.character(g)[match(d$sample,samples)],levels=pops);d}
 popscale <- function()list(scale_color_manual(values=popcols,drop=FALSE),scale_shape_manual(values=popshape,drop=FALSE))
 palout <- list()
 # Family-only population means, with a fixed shared categorical palette in both assays.
 txall<-relative(t(rowsum(C,tax$Family,reorder=TRUE)));ix<-prot$superKingdom=='Bacteria';pxall<-relative(t(rowsum(V[ix,,drop=FALSE],prot$Family[ix],reorder=TRUE)))
 taxa<-sort(union(colnames(txall),colnames(pxall)));a<-b<-setNames(rep(0,length(taxa)),taxa);a[colnames(txall)]<-colMeans(txall);b[colnames(pxall)]<-colMeans(pxall)
 ranked<-names(sort(pmax(a,b),decreasing=TRUE));tops<-head(setdiff(ranked,'unclassified'),12);display<-c(tops,'Other','unclassified')
 tc<-setNames(c('#332288','#88CCEE','#44AA99','#117733','#999933','#DDCC77','#CC6677','#882255','#AA4499','#661100','#6699CC','#888888')[seq_along(tops)],tops);tc<-c(tc,Other='#BDBDBD',unclassified='#EEEEEE')
 labels<-setNames(display,display);labels[tops]<-ifelse(tops%in%ca$Family,paste0(tops,' *'),tops);palout<-list(Family=data.frame(rank='Family',taxon=names(tc),hex=unname(tc)))
 for(cy in c('with_Cyanobacteria','without_Cyanobacteria')){
  if(cy=='with_Cyanobacteria'){tx<-txall;px<-pxall} else {keep<-tax$Phylum!='Cyanobacteria';tx<-relative(t(rowsum(C[keep,,drop=FALSE],tax$Family[keep],reorder=TRUE)));keep<-prot$partner=='Noncyano_bacteria';px<-relative(t(rowsum(V[keep,,drop=FALSE],prot$Family[keep],reorder=TRUE)))}
  plotlist<-list()
  for(kind in c('16S','Protein')){
   x<-if(kind=='16S')tx else px;z<-matrix(0,nrow(x),length(display),dimnames=list(rownames(x),display));present<-intersect(tops,colnames(x));z[,present]<-x[,present,drop=FALSE];rest<-setdiff(colnames(x),c(tops,'unclassified'));z[,'Other']<-rowSums(x[,rest,drop=FALSE]);if('unclassified'%in%colnames(x))z[,'unclassified']<-x[,'unclassified']
   stopifnot(max(abs(rowSums(z)-1))<1e-8,all(display%in%names(tc)),!anyNA(tc[display]))
   d<-long(z);d$feature<-factor(d$feature,levels=display);dm<-aggregate(value~population+feature,d,mean)
   plotlist[[kind]]<-ggplot(dm,aes(population,value*100,fill=feature))+geom_col(width=.8,color='white',linewidth=.1,show.legend=TRUE)+scale_fill_manual(values=tc,limits=display,breaks=display,labels=function(x)unname(labels[x]),name='Family',drop=FALSE)+scale_y_continuous(expand=expansion(mult=c(0,.01)),limits=c(0,100))+labs(x=NULL,y=if(kind=='16S')'Bacterial reads (%)' else 'Bacterial NSAF (%)',title=if(kind=='16S')'A  16S families' else 'B  Bacterial protein families')+theme(axis.text.x=element_text(angle=50,hjust=1,size=9),plot.title=element_text(size=10),legend.position='right',legend.text=element_text(size=8),legend.title=element_text(size=9),legend.key.height=grid::unit(4,'mm'),legend.key.width=grid::unit(4.5,'mm'),legend.spacing.y=grid::unit(0,'mm'),legend.margin=margin(0,0,0,0))+guides(fill=guide_legend(ncol=1,byrow=TRUE,override.aes=list(alpha=1)))
  }
  savefig(combine_plots(plotlist,nrow=1,guides='collect',legend_position='right'),paste0('Taxonomy_Family_',cy),paste('Population means of three thalli each; all cleaned bacterial reads or bacterial NSAF are included before display aggregation.',gsub('_',' ',cy),'. The same top 12 categories, colors, and legend are retained in both assays and exclusion sensitivity; other named families form Other. Various denotes unresolved/multiple family assignments, not one biological family. Unclassified is retained in the denominator. Asterisk denotes a family containing at least one strict Kilimanjaro-core ASV; it does not assign proteins to individual ASVs. Read and NSAF percentages measure different quantities.'),width=10.5,height=4.8)
 }
 write.csv(do.call(rbind,palout),file.path(out,'taxonomy_color_palette.csv'),row.names=FALSE)
 write.csv(data.frame(population=pops,hex=popcols,shape=popshape),file.path(out,'population_color_palette.csv'),row.names=FALSE)
 # Full strict-core family inventory (including unresolved assignments).
 coretab<-aggregate(cbind(Kilimanjaro=rep(1,nrow(ca)),Styria_all3=as.integer(ca$Styria_detections==3)),list(Family=ca$Family),sum)
 coretab<-coretab[order(coretab$Kilimanjaro),];cm<-as.matrix(coretab[,-1]);rownames(cm)<-coretab$Family;ct<-as.data.frame(as.table(cm),stringsAsFactors=FALSE);names(ct)<-c('family','criterion','count');ct$family<-factor(ct$family,levels=coretab$Family)
 p<-ggplot(ct,aes(criterion,family,fill=count))+geom_tile(color='white',linewidth=.5)+geom_text(aes(label=count,color=count>max(count)*.55),size=3.5)+scale_color_manual(values=c('FALSE'='black','TRUE'='white'),guide='none')+scale_fill_gradient(low='#EFF3FF',high='#08519C',name='Core ASVs')+scale_x_discrete(labels=c(Kilimanjaro='Strict Kilimanjaro core\n18/18 thalli',Styria_all3='Recovered in Styria\n3/3 thalli'))+labs(x=NULL,y=NULL,title='Strict-core ASV recovery',subtitle=sprintf('%d strict Kilimanjaro-core ASVs; %d recovered in all three Styria thalli',nrow(ca),sum(ca$Styria_detections==3)))+theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position='right')
 savefig(p,'Core_family_recovery','The strict Kilimanjaro core comprises 118 ASVs detected with at least one read in all 18 Kilimanjaro thalli. Of these, 92 were recovered in all three Styria thalli and hence detected in all 21 sampled thalli. Styria was excluded from initial core identification. The two columns show family-level counts of these ASV sets, not family-level core definitions or counts of strains. Recovery in at least one Styria thallus is reported separately in results tables and is not displayed as a strict shared core.',height=max(5.5,nrow(coretab)*.23+1.5))
 # Display the 25 most abundant strict-core ASVs; complete 118-ASV table is retained.
 coreids<-ca$ASV[order(ca$mean_Kilimanjaro_fraction,decreasing=TRUE)];coreids<-head(coreids,25)
 cm<-t(C[coreids,,drop=FALSE]>0)*1;pd<-long(cm);labsasv<-setNames(paste0(coreids,' | ',tax[coreids,'Family']),coreids);pd$feature<-factor(pd$feature,levels=rev(coreids));pd$sample<-factor(pd$sample,levels=samples)
 p<-ggplot(pd,aes(sample,feature,fill=factor(value)))+geom_tile(color='white',linewidth=.3)+scale_fill_manual(values=c('0'='#EEEEEE','1'='#0072B2'),name='Detection',labels=c('Not detected','Detected'),drop=FALSE)+scale_y_discrete(labels=labsasv)+facet_grid(~population,scales='free_x',space='free_x')+labs(x='Individual thallus',y=NULL,title='Detection of abundant strict-core ASVs',subtitle='Top 25 by mean Kilimanjaro read fraction; full inventory retained')+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5),axis.text.y=element_text(size=8),strip.text=element_text(size=8),axis.line=element_blank(),axis.ticks=element_blank())
 savefig(p,'Core_ASV_detection','Illustrative subset: 25 strict-core ASVs selected by mean Kilimanjaro relative abundance. Presence requires at least one read; absence is non-detection at the observed sequencing depth. All 118 strict-core ASVs and Styria recovery counts remain in core_ASVs.csv. Family labels do not establish strain identity.',height=8.5)

 # Family-level alpha diversity, with/without Cyanobacteria: one value per thallus.
 # Unresolved bins are excluded from named-family diversity and their read share is reported.
 set.seed(20261101);nr<-200L;family_counts<-list();family_exclusions<-list()
 for(co in c('All_bacteria','Noncyano_bacteria')){
  keepasv<-if(co=='All_bacteria')rep(TRUE,nrow(C)) else tax$Phylum!='Cyanobacteria'
  famlabel<-trimws(as.character(tax$Family));resolved<-!tolower(famlabel)%in%c('unclassified','various','unknown','na','','unassigned')
  z<-t(rowsum(C[keepasv & resolved,,drop=FALSE],famlabel[keepasv & resolved],reorder=TRUE));family_counts[[co]]<-z
  family_exclusions[[co]]<-data.frame(compartment=co,sample=samples,total_compartment_reads=colSums(C[keepasv,,drop=FALSE]),named_family_reads=rowSums(z),excluded_unresolved_reads=colSums(C[keepasv & !resolved,,drop=FALSE]))
 }
 depthalpha<-min(vapply(family_counts,function(x)min(rowSums(x)),numeric(1)));if(depthalpha<1)stop('No named-family reads in one or more alpha-diversity profiles.')
 alpha_one<-function(x)data.frame(sample=rownames(x),Observed_families=rowSums(x>0),Shannon=vegan::diversity(x,index='shannon'),Inverse_Simpson=vegan::diversity(x,index='invsimpson'))
 rawalpha<-list();drawalpha<-list()
 for(co in names(family_counts)){
  z<-family_counts[[co]];a<-alpha_one(z);a$compartment<-co;rawalpha[[co]]<-a
  for(b in seq_len(nr)){a<-alpha_one(vegan::rrarefy(z,sample=depthalpha));a$compartment<-co;a$draw<-b;drawalpha[[length(drawalpha)+1L]]<-a}
 }
 rawalpha<-do.call(rbind,rawalpha);ad<-do.call(rbind,drawalpha)
 al<-reshape(ad,direction='long',varying=c('Observed_families','Shannon','Inverse_Simpson'),v.names='value',timevar='metric',times=c('Observed_families','Shannon','Inverse_Simpson'),idvar=c('sample','compartment','draw'))
 am<-aggregate(value~sample+compartment+metric,al,mean);am$population<-factor(as.character(g)[match(am$sample,samples)],levels=pops)
 am$Cyanobacteria<-factor(ifelse(am$compartment=='All_bacteria','Including Cyanobacteria','Excluding Cyanobacteria'),levels=c('Including Cyanobacteria','Excluding Cyanobacteria'))
 am$metric<-factor(am$metric,levels=c('Observed_families','Shannon','Inverse_Simpson'),labels=c('Observed named families','Shannon diversity','Inverse Simpson diversity'))
 write.csv(rawalpha,file.path(out,'alpha_family_unrarefied.csv'),row.names=FALSE);write.csv(ad,file.path(out,'alpha_family_rarefaction_draws.csv'),row.names=FALSE);write.csv(am,file.path(out,'alpha_family_thallus_means.csv'),row.names=FALSE);write.csv(do.call(rbind,family_exclusions),file.path(out,'alpha_family_read_exclusions.csv'),row.names=FALSE)
 alpha_plot<-function(d)ggplot(d,aes(population,value,fill=population))+geom_boxplot(width=.65,outlier.shape=NA,alpha=.65,linewidth=.45)+geom_point(aes(shape=population),position=position_jitter(width=.10,height=0,seed=1),size=2.4,color='black')+scale_fill_manual(values=popcols)+scale_shape_manual(values=popshape)+labs(x=NULL,y='Family-level alpha diversity',subtitle=paste('Three thalli per population;',nr,'rarefactions per thallus at',depthalpha,'named-family reads'))+theme(axis.text.x=element_text(angle=50,hjust=1,size=9),legend.position='none')
 captionalpha<-paste('Family-level 16S alpha diversity including and excluding Cyanobacteria. Unresolved family bins are excluded and their read counts exported. Each dot is one thallus mean across',nr,'rarefactions at a common depth of',depthalpha,'reads across both compartments. Each box contains three biological replicates, not 200 rarefaction draws: center line median, hinges 25th/75th percentiles, whiskers within 1.5 times IQR. No significance stars or population tests are added. Raw metrics and all rarefaction draws are retained.')
 p<-alpha_plot(am)+facet_grid(metric~Cyanobacteria,scales='free_y')+labs(title='Bacterial family-level alpha diversity')
 savefig(p,'Alpha_family_boxplots_with_without_Cyanobacteria',captionalpha,width=7.2,height=8.2)
 
 qc<-data.frame(population=g,sample=samples,reads=colSums(C));p<-ggplot(qc,aes(population,reads,color=population,shape=population))+geom_point(size=3,position=position_jitter(width=.08,height=0,seed=1))+popscale()+labs(x=NULL,y='Clean bacterial reads',title='Sequencing depth by population')+theme(axis.text.x=element_text(angle=45,hjust=1),legend.position='none');savefig(p,'QC_sequencing_depth','Each point is one thallus. Counts contextualize detection-based core and alpha-diversity analyses.',height=4.5)
 # Three panels per row, for OG/subRole and with/without Cyanobacteria.
 nmcached<-list();nmstress<-list()
 nmpanel<-function(x,key,title){
  if(is.null(nmcached[[key]])){
   D<-vegan::vegdist(x,method='bray');set.seed(20261101);nm<-vegan::metaMDS(D,k=2,trymax=100,trace=FALSE,autotransform=FALSE)
   xy<-as.data.frame(vegan::scores(nm,display='sites'));names(xy)<-c('NMDS1','NMDS2');xy$population<-g;xy$sample<-samples
   nmcached[[key]]<<-xy;nmstress[[key]]<<-nm$stress;write.csv(xy,file.path(out,paste0('NMDS_',key,'.csv')),row.names=FALSE)
  }
  xy<-nmcached[[key]];cx<-mean(range(xy$NMDS1));cyy<-mean(range(xy$NMDS2));half<-max(diff(range(xy$NMDS1)),diff(range(xy$NMDS2)))*.60;if(half<=0)half<-.1
  ggplot(xy,aes(NMDS1,NMDS2,fill=population))+geom_point(shape=21,size=3.6,color='black',stroke=.35,alpha=.9)+scale_fill_manual(values=popcols,drop=FALSE,name='Population')+coord_equal(xlim=c(cx-half,cx+half),ylim=c(cyy-half,cyy+half),expand=FALSE)+labs(x='NMDS1',y='NMDS2',title=title,subtitle=sprintf('Bray-Curtis; stress = %.3f',nmstress[[key]]))+theme(plot.title=element_text(size=10),axis.title=element_text(size=10),axis.title.x=element_text(color=if(startsWith(title,'A'))'black' else 'transparent'),axis.text=element_text(size=8),aspect.ratio=1,panel.border=element_rect(color='black',fill=NA,linewidth=.6),axis.line=element_blank())+guides(fill=guide_legend(nrow=1,override.aes=list(shape=21,size=3.6,color='black',alpha=1)))

 }
 for(cy in c('with_Cyanobacteria','without_Cyanobacteria'))for(le in c('og','subRole')){
  co<-if(cy=='with_Cyanobacteria')'All_bacteria' else 'Noncyano_bacteria'
  tx<-profiles[[paste('16S',co,'Family',sep='__')]];bx<-profiles[[paste('Function',co,le,sep='__')]]
  wx<-if(cy=='with_Cyanobacteria')profiles[[paste('Function','Whole',le,sep='__')]] else relative(t(rowsum(V[prot$partner!='Cyanobiont',,drop=FALSE],prot[[le]][prot$partner!='Cyanobiont'],reorder=TRUE)))
  ps<-list(nmpanel(tx,paste('16S',co,'Family',sep='__'),'A  16S families'),nmpanel(bx,paste('Bacteria',co,le,sep='__'),paste('B  Bacterial metaproteome:',le)),nmpanel(wx,paste('Whole',cy,le,sep='__'),paste('C  Whole metaproteome:',le)))
  savefig(combine_plots(ps,nrow=1,guides='collect'),paste('NMDS_ABC',le,cy,sep='_'),paste('A: family-level 16S; B: bacterial functional profiles; C: whole-metaproteome functional profiles.',gsub('_',' ',cy),'. Excluding Cyanobacteria removes their reads from A, cyanobacterial proteins from B, and cyanobacterial proteins only from C before renormalization; fungal and algal partners remain in C. One filled circle per thallus; population is encoded by fill color only, with no connecting segments. One shared population legend and NMDS1 title on the leftmost panel only. Panels have independent NMDS configurations and axes; compare within-panel structure, not coordinates between panels. No ellipses are fitted with n=3 per population. NMDS geometry is not used in centroid-distance inference.'),width=10.5,height=4.4)
 }
 # Exploratory taxon-correlation plots omitted from the manuscript figure set.
 # Mantel matrix summary: reported coefficients and adjusted P values from actual tests.
 md<-mm;md$label<-sprintf('%.2f%s',md$r,ifelse(md$q<.05,' *',''));md$scope<-factor(md$scope,levels=unique(md$scope));md$function_level<-factor(md$function_level,levels=c('og','subRole','mainRole'))
 p<-ggplot(md,aes(function_level,taxonomy_level,fill=r))+geom_tile(color='white',linewidth=.5)+geom_text(aes(label=label),size=3.5)+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-1,1),name='Mantel r')+facet_grid(scope~compartment,labeller=label_wrap_gen(20))+labs(x='Functional resolution',y='Taxonomic resolution',title='D  Taxonomic–functional correspondence')+theme(axis.line=element_blank(),axis.ticks=element_blank())
 savefig(p,'Mantel_test_summary','Cells show two-sided Pearson Mantel r between Bray–Curtis dissimilarity matrices. Asterisk: BH-adjusted q < 0.05 across 18 comparisons within each scope. Scope and bacterial compartment are explicit; tests are separate from the feature-correlation plots.',height=5.8)
 for(co in c('All_bacteria','Noncyano_bacteria'))for(sc in c('all_seven','Kilimanjaro')){
  ids<-if(sc=='all_seven')seq_along(samples) else which(g!='Styria');x<-as.vector(vegan::vegdist(profiles[[paste('16S',co,'Family',sep='__')]][ids,,drop=FALSE]));y<-as.vector(vegan::vegdist(profiles[[paste('Function',co,'og',sep='__')]][ids,,drop=FALSE]));d<-data.frame(taxonomy=x,functional=y);test<-mm[mm$scope==sc & mm$compartment==co & mm$taxonomy_level=='Family' & mm$function_level=='og',]
  if(nrow(test)!=1)stop('Mantel plotting key does not identify one test: ',co,' ',sc)
  p<-ggplot(d,aes(taxonomy,functional))+geom_point(shape=16,size=1.5,alpha=.5,color='#3B528B')+labs(x='16S family Bray–Curtis dissimilarity',y='Bacterial OG Bray–Curtis dissimilarity',title=gsub('_',' ',co),subtitle=sprintf('%s | Mantel r = %.3f; P = %.4g; q = %.4g',sc,test$r,test$p,test$q))
  savefig(p,paste('Mantel_scatter',co,sc,sep='_'),'Each point is a pair of thalli. Pairwise distances are not independent observations; inference comes from the matrix-permutation Mantel test shown, not from an ordinary regression. No independent-point confidence band or regression P value is drawn.',height=4.8)
 }
 # Centroid matrices at both functional resolutions and both bacterial definitions.
 centroid_cache<-list()
 for(cy in c('with_Cyanobacteria','without_Cyanobacteria'))for(le in c('og','subRole')){
  co<-if(cy=='with_Cyanobacteria')'All_bacteria' else 'Noncyano_bacteria'
  DS<-list(cent[[paste('16S',co,'Family',sep='__')]],cent[[paste('Function',co,le,sep='__')]])
  if(cy=='with_Cyanobacteria')DS[[3]]<-cent[[paste('Function','Whole',le,sep='__')]] else {
   wx<-relative(t(rowsum(V[prot$partner!='Cyanobiont',,drop=FALSE],prot[[le]][prot$partner!='Cyanobiont'],reorder=TRUE)))
   pc<-vegan::wcmdscale(vegan::vegdist(wx,method='bray'),k=nrow(wx)-1,eig=TRUE,add='lingoes');Z<-pc$points
   CT<-t(vapply(pops,function(pop)colMeans(Z[g==pop,,drop=FALSE]),numeric(ncol(Z))));DS[[3]]<-as.matrix(dist(CT));dimnames(DS[[3]])<-list(pops,pops)
  }
  vmax<-max(vapply(DS,max,numeric(1)));ps<-list()
  for(i in 1:3){d<-as.data.frame(as.table(DS[[i]]),stringsAsFactors=FALSE);names(d)<-c('from','to','distance');d$from<-factor(d$from,levels=pops);d$to<-factor(d$to,levels=rev(pops));ps[[i]]<-ggplot(d,aes(from,to,fill=distance))+geom_tile(color='white',linewidth=.2)+scale_fill_viridis_c(limits=c(0,vmax),name='Centroid distance')+coord_equal()+labs(x=NULL,y=NULL,title=c('A  16S families',paste('B  Bacterial',le),paste('C  Whole-metaproteome',le))[i])+theme(axis.text.x=element_text(angle=45,hjust=1),axis.text=element_text(size=8),axis.line=element_blank(),axis.ticks=element_blank(),plot.title=element_text(size=10),legend.position='right',legend.direction='vertical',legend.title=element_text(angle=90))+guides(fill=guide_colorbar(direction='vertical',barheight=grid::unit(34,'mm'),barwidth=grid::unit(3.5,'mm')))}
  savefig(combine_plots(ps,nrow=1,guides='collect',legend_position='right'),paste('Centroid_ABC',le,cy,sep='_'),paste('A: 16S families; B: bacterial functions; C: whole-metaproteome functions.',gsub('_',' ',cy),'. Bray-Curtis dissimilarities were embedded in all positive PCoA dimensions after Lingoes correction; Euclidean distances between population means in that space are shown. One common color scale within this figure. Each population contains three thalli. Exclusion removes only cyanobacterial reads/proteins before normalization. These are distances, not correlations or NMDS distances. At OG resolution including Cyanobacteria, each Kilimanjaro centroid is farther from Styria than from each other Kilimanjaro centroid. The separate exact reference-relabelling test has seven choices and minimum attainable P = 1/7 = 0.143; conventional outlier significance is unavailable.'),width=10.5,height=4.1)
 }
 # Direction labels sit in a colored strip immediately above each volcano panel.
 volcano_direction_strip<-function(p,left_label,right_label){
  gt<-ggplot2::ggplotGrob(p);ii<-which(gt$layout$name=='panel');if(length(ii)!=1L)stop('Expected one volcano panel.')
  panel<-gt$layout[ii,];gt<-gtable::gtable_add_rows(gt,grid::unit(5.5,'mm'),pos=panel$t-1L)
  strip<-grid::grobTree(grid::rectGrob(x=.25,width=.5,gp=grid::gpar(fill='#0072B2',col=NA)),grid::rectGrob(x=.75,width=.5,gp=grid::gpar(fill='#D55E00',col=NA)),grid::textGrob(left_label,x=.25,gp=grid::gpar(col='white',fontsize=8,fontface='bold')),grid::textGrob(right_label,x=.75,gp=grid::gpar(col='white',fontsize=8,fontface='bold')))
  gtable::gtable_add_grob(gt,strip,t=panel$t,l=panel$l,r=panel$r,clip='off',name='direction-strip')
 }
 # Select an extreme population by a single explicit rule, before feature-level display.
 primaryD<-cent[['Function__All_bacteria__og']];candidate<-pops[pops!='Styria'];selected<-candidate[which.max(primaryD[candidate,'Styria'])]
 selection<-data.frame(population=candidate,bacterial_OG_distance_to_Styria=primaryD[candidate,'Styria'],selected=candidate==selected,selection_rule='Maximum full-space all-bacterial OG centroid distance to Styria; same population used for every volcano panel')
 write.csv(selection,file.path(out,'volcano_extreme_population_selection.csv'),row.names=FALSE)
 cat('\nExploratory volcano population selected by bacterial OG centroid distance: ',selected,'\n',sep='')
 da16<-read.csv(file.path(out,'DA','taxonomy_pairwise_DESeq2.csv'),check.names=FALSE);daf<-read.csv(file.path(out,'DA','functional_pairwise_limma.csv'),check.names=FALSE)
 # Additional whole-metaproteome OG model excluding cyanobacterial protein groups.
 # Same detected-log model, detection eligibility and 21-pair correction as native inference.
 whole_no_cy_limma <- function(){
  needed<-c('limma','statmod');if(!all(vapply(needed,requireNamespace,logical(1),quietly=TRUE)))stop('Install limma and statmod for the whole-metaproteome exclusion volcano panel.')
  mask<-prot$partner!='Cyanobiont';v<-V[mask,,drop=FALSE];M<-rowsum(v,prot$og[mask],reorder=FALSE);M<-sweep(M,2,colSums(v),'/');M<-M[!rownames(M)%in%c('unclassified','various','Function unknown','poorly characterized'),,drop=FALSE]
  Y<-log2(M);Y[!is.finite(Y)]<-NA;nd<-sapply(pops,function(pop)rowSums(M[,g==pop,drop=FALSE]>0));rownames(nd)<-rownames(M);df<-rowSums(is.finite(Y))-rowSums(nd>0);keep<-rowSums(is.finite(Y))>=6 & df>=3
  if(sum(keep)<2)stop('Too few eligible whole-metaproteome OGs after excluding Cyanobacteria.')
  design<-model.matrix(~0+g);colnames(design)<-pops;fit<-limma::lmFit(Y[keep,,drop=FALSE],design);ans<-list()
  for(pair in combn(pops,2,simplify=FALSE)){
   ids<-rownames(M)[keep & nd[,pair[1]]>=2 & nd[,pair[2]]>=2];if(length(ids)<2)next
   ct<-matrix(0,length(pops),1,dimnames=list(pops,NULL));ct[pair[1],1]<- -1;ct[pair[2],1]<-1;ff<-limma::contrasts.fit(fit[match(ids,rownames(fit$coefficients)),],ct);n<-nrow(ff$coefficients)
   ff<-tryCatch(limma::eBayes(ff,trend=n>=30,robust=n>=20),error=function(e){warning('Whole-exclusion empirical-Bayes fallback: ',conditionMessage(e));limma::eBayes(ff,trend=FALSE,robust=FALSE)})
   tt<-limma::topTable(ff,number=Inf,sort.by='none')
   ans[[length(ans)+1L]]<-data.frame(compartment='Whole_without_Cyanobacteria',level='og',popA=pair[1],popB=pair[2],feature=rownames(tt),log2FoldChange=tt$logFC,pvalue=tt$P.Value,n_detected_A=nd[rownames(tt),pair[1]],n_detected_B=nd[rownames(tt),pair[2]])
  }
  z<-do.call(rbind,ans);if(is.null(z)||!nrow(z))stop('No eligible whole-metaproteome exclusion contrasts.')
  z$q_across_21_contrasts<-p.adjust(z$pvalue,'BH');write.csv(z,file.path(out,'DA','whole_without_Cyanobacteria_OG_pairwise_limma.csv'),row.names=FALSE);z
 }
 whole_exclusion<-whole_no_cy_limma()
 for(cy in c('with_Cyanobacteria','without_Cyanobacteria')){
  co<-if(cy=='with_Cyanobacteria')'All_bacteria' else 'Noncyano_bacteria';ps<-list();ex<-list()
  for(kind in c('16S','Bacterial_OG','Whole_OG')){
   z<-if(kind=='16S')da16[da16$compartment==co & da16$level=='Family' & da16$popA==selected & da16$popB=='Styria',] else if(kind=='Bacterial_OG')daf[daf$compartment==co & daf$level=='og' & daf$popA==selected & daf$popB=='Styria',] else if(cy=='with_Cyanobacteria')daf[daf$compartment=='Whole' & daf$level=='og' & daf$popA==selected & daf$popB=='Styria',] else whole_exclusion[whole_exclusion$popA==selected & whole_exclusion$popB=='Styria',]
   z<-z[is.finite(z$log2FoldChange)&is.finite(z$q_across_21_contrasts),];if(!nrow(z))stop('No finite native contrast for selected population: ',selected,' ',kind,' ',co)
   z$mlogq<--log10(pmax(z$q_across_21_contrasts,1e-300));z$status<-ifelse(z$q_across_21_contrasts>=.05,'Not significant',ifelse(z$log2FoldChange<0,paste('Higher in',selected),'Higher in Styria'));z$status<-factor(z$status,levels=c('Not significant',paste('Higher in',selected),'Higher in Styria'));cols<-setNames(c('#BDBDBD','#0072B2','#D55E00'),levels(z$status));sig<-z[z$q_across_21_contrasts<.05,];lab<-head(sig[order(sig$q_across_21_contrasts),],6);ex[[kind]]<-z
   vp<-ggplot(z,aes(log2FoldChange,mlogq,color=status))+geom_hline(yintercept=-log10(.05),linetype=2,linewidth=.4,color='grey45')+geom_vline(xintercept=0,linewidth=.3,color='grey60')+geom_point(size=2.5,alpha=.7)+scale_color_manual(values=cols,drop=FALSE,guide='none')+labs(x=paste0('Model log2(Styria / ',selected,')'),y=if(kind=='16S')'-log10(BH-adjusted P)' else NULL,title=switch(kind,'16S'='A  16S families',Bacterial_OG='B  Bacterial OGs',Whole_OG='C  Whole-metaproteome OGs'),subtitle=sprintf('%s vs Styria: %d / %d significant',selected,nrow(sig),nrow(z)))+theme(legend.position='none',plot.title=element_text(size=10),plot.subtitle=element_text(size=8),axis.title=element_text(size=10),axis.text=element_text(size=8),aspect.ratio=1,panel.border=element_rect(color='black',fill=NA,linewidth=.6),axis.line=element_blank())
   ps[[kind]]<-volcano_direction_strip(vp,paste('Higher in',selected),'Higher in Styria')

   write.csv(z,file.path(out,paste('volcano_selected',kind,cy,'csv',sep='.')),row.names=FALSE)
  }
  savefig(combine_plots(ps,nrow=1),paste0('Volcano_extreme_',cy),paste('Styria versus',selected,'selected as the farthest population in all-bacterial OG full-space centroid distance. Selection is exploratory and data-dependent; it does not provide an independent confirmatory comparison. The same population is retained in the exclusion sensitivity. A: native DESeq2 family-count contrasts; B: bacterial OGs; C: whole-metaproteome OGs, both using limma on detected log-relative NSAF. The whole-metaproteome exclusion panel is fitted separately after cyanobacterial protein removal using the same seven-population model, eligibility rules and correction across all 21 pairs. All eligible finite-q features are shown without feature labels, with point size 2.5 and alpha 0.7; the shared y-axis title is printed on panel A only. Direction banners sit above each panel. Gray nonsignificant, blue higher in the selected Kilimanjaro population, orange higher in Styria. Dashed line q = 0.05; no additional fold-change cutoff. BH correction retains all eligible features and all 21 population pairs within compartment/resolution, not only the displayed contrast. Relative protein representation is not absolute activity.'),width=11.5,height=4.8)
 }
 # Descriptive ASV and family differences: popB minus popA, proportions in percentage points.
 allpairs<-combn(pops,2,simplify=FALSE);relasv<-sweep(C,2,colSums(C),'/');asvmeans<-sapply(pops,function(pop)rowMeans(relasv[,g==pop,drop=FALSE]));rownames(asvmeans)<-rownames(C)
 asvsummary<-data.frame(ASV=rownames(C),Family=tax$Family,asvmeans,check.names=FALSE);write.csv(asvsummary,file.path(out,'ASV_population_mean_read_fractions.csv'),row.names=FALSE)
 asvdif<-do.call(rbind,lapply(allpairs,function(pair)data.frame(ASV=rownames(C),Family=tax$Family,popA=pair[1],popB=pair[2],mean_fraction_A=asvmeans[,pair[1]],mean_fraction_B=asvmeans[,pair[2]],difference_percentage_points=100*(asvmeans[,pair[2]]-asvmeans[,pair[1]]))))
 write.csv(asvdif,file.path(out,'ASV_pairwise_mean_differences_descriptive.csv'),row.names=FALSE)
 famrel<-relative(t(rowsum(C,tax$Family,reorder=TRUE)));fm<-sapply(pops,function(pop)colMeans(famrel[g==pop,,drop=FALSE]));rownames(fm)<-colnames(famrel)
 fd<-do.call(rbind,lapply(allpairs,function(pair)data.frame(feature=rownames(fm),popA=pair[1],popB=pair[2],mean_fraction_A=fm[,pair[1]],mean_fraction_B=fm[,pair[2]],difference_percentage_points=100*(fm[,pair[2]]-fm[,pair[1]]))))
 native<-da16[da16$compartment=='All_bacteria' & da16$level=='Family',c('feature','popA','popB','log2FoldChange','pvalue','q_across_21_contrasts')]
 fd<-merge(fd,native,by=c('feature','popA','popB'),all.x=TRUE,sort=FALSE);write.csv(fd,file.path(out,'Family_pairwise_mean_differences_and_native_DESeq2.csv'),row.names=FALSE)
 write.csv(fd[fd$feature=='Sphingomonadaceae',],file.path(out,'Sphingomonadaceae_16S_pairwise_summary.csv'),row.names=FALSE)
 # Descriptive mean ASV differences carry no invented ASV-level P values.
 # Broad partner functions, detailed OGs, and family contributions.
 for(le in c('subRole','og')){
  panels<-list()
  for(co in names(partner_names)){
   z<-rankings[rankings$partner==co & rankings$level==le & rankings$informative_annotation,];z<-head(z[order(z$rank_informative),],if(le=='og')8 else 8)
   z$label<-if(le=='og')paste(z$feature,substr(z$description,1,100),sep=': ') else z$feature;z$label<-factor(wrap(z$label,36),levels=rev(wrap(z$label,36)))
   panels[[co]]<-ggplot(z,aes(mean_fraction*100,label))+geom_col(fill=partnercols[co],width=.7)+labs(x='Mean within-partner NSAF (%)',y=NULL,title=paste(LETTERS[match(co,names(partner_names))],partner_names[co]))+scale_x_continuous(expand=expansion(mult=c(0,.06)))+theme(axis.text.y=element_text(size=9),plot.title=element_text(size=11))
  }
  savefig(combine_plots(panels,ncol=2),paste0('Partner_functions_2x2_',le),'Panels: cyanobiont, mycobiont, non-cyanobacterial bacteria, phycobiont. Top eight informative functions per partner; unknown bins remain in denominators. Means equally weight 21 thalli. Relative protein representation is not metabolic flux; descriptions are database annotations. Each panel has its own abundance axis.',width=11,height=if(le=='og')11 else 8.6)
 }
 fc<-do.call(rbind,coverage);eligible<-ft[ft$level=='og' & ft$eligible_all_seven,]$family;fc<-fc[fc$family%in%eligible,];fc$population<-factor(fc$population,levels=pops)
 p<-ggplot(fc,aes(population,all_bacterial_fraction*100,color=population,shape=population))+geom_point(position=position_jitter(width=.12,height=0,seed=1),size=2.2)+popscale()+facet_wrap(~family,scales='free_y',ncol=2)+labs(x=NULL,y='Family share of bacterial NSAF (%)',title='A  Shared-family protein contributions')+theme(axis.text.x=element_text(angle=45,hjust=1,size=8),legend.position='none')
 savefig(p,'Shared_family_protein_contributions','Each point is one thallus. Families are linked by taxonomic annotation, not by exact ASV-to-protein assignment. Panels use different y scales; detection eligibility and native contribution tests are reported in complete tables.',height=8.8)
 z<-ft[ft$eligible_all_seven,];z$label<-sprintf('q = %.3g',z$q);p<-ggplot(z,aes(R2,reorder(family,R2),color=level,shape=level))+geom_point(size=2.7,position=position_dodge(width=.5))+geom_text(aes(label=label),position=position_dodge(width=.5),hjust=-.15,size=2.7,show.legend=FALSE)+scale_color_manual(values=c(og='#0072B2',subRole='#D55E00'))+scale_x_continuous(limits=c(0,max(z$R2)+.22))+labs(x='Within-family population PERMANOVA R-squared',y=NULL,color='Resolution',shape='Resolution',title='B  Within-family functional differentiation')
 savefig(p,'Shared_family_functional_structure','BH-adjusted q covers the 14 eligible family-by-resolution tests. Each family requires at least two positive protein profiles per population. Nonsignificance is not equivalence; family-level annotations do not establish functions of individual shared ASVs.',height=5.6)
 # Core-definition and sequencing-depth sensitivities remain explicitly available.
 thpath<-file.path(out,'core_thresholds.csv');rapath<-file.path(out,'core_rarefaction.csv')
 if(!all(file.exists(c(thpath,rapath))))stop('Core sensitivity tables missing; run RUN_MODE = full to regenerate.')
 th<-read.csv(thpath,check.names=FALSE);ra<-read.csv(rapath,check.names=FALSE)
 th$definition<-factor(th$definition,levels=c('one_read','two_reads','relative_0.01pct'),labels=c('At least 1 read','At least 2 reads','At least 0.01% reads'))
 pa<-ggplot(th,aes(factor(min_Kilimanjaro_detections),core_ASVs,fill=definition))+geom_col(position='dodge',width=.75)+facet_wrap(~compartment,nrow=1)+scale_fill_manual(values=c('#0072B2','#D55E00','#009E73'),name='Detection threshold')+labs(x='Minimum Kilimanjaro detections (18 thalli)',y='Core ASVs',title='A  Sensitivity to core definition')+theme(legend.position='bottom',legend.text=element_text(size=8))
 rb<-reshape(ra[,c('draw','core_ASVs','Styria_all3')],direction='long',varying=c('core_ASVs','Styria_all3'),v.names='ASVs',timevar='criterion',times=c('Kilimanjaro strict core','Recovered in Styria 3/3'),idvar='draw')
 pb<-ggplot(rb,aes(criterion,ASVs,fill=criterion))+geom_boxplot(outlier.shape=NA,width=.65)+scale_fill_manual(values=c('#0072B2','#D55E00','#009E73'))+labs(x=NULL,y='ASVs per rarefaction draw',title='B  Rarefaction sensitivity',subtitle=paste(nrow(ra),'draws;',unique(ra$depth),'reads per thallus'))+theme(legend.position='none',axis.text.x=element_text(angle=20,hjust=1))
 savefig(combine_plots(list(pa,pb),ncol=1),'Core_definition_sensitivity','A: prevalence and detection-threshold sensitivities for bacterial compartments. B: variation among 200 rarefaction draws; these draws are computational subsamples, not biological replicates. The main core remains >=1 read in all 18 Kilimanjaro thalli. Styria is used only to evaluate recovery.',height=8)
 # Numbering is explicit and independent of loop or execution order.
 numbered <- list()
 add_numbered <- function(id,title,keys,plot=NULL,width=NULL,height=NULL){
  stopifnot(all(keys %in% names(plotcache)))
  entries<-plotcache[keys]
  if(is.null(plot))plot<-entries[[1]]$plot
  if(is.null(width))width<-entries[[1]]$width
  if(is.null(height))height<-entries[[1]]$height
  numbered[[id]] <<- list(title=title,plot=plot,caption=paste(vapply(entries,function(x)x$caption,character(1)),collapse='\n\n'),width=width,height=height,source=paste(keys,collapse='; '))
 }
 getplot<-function(key)plotcache[[key]]$plot
 add_numbered('Figure1','Bacterial family composition: reads and proteins','Taxonomy_Family_with_Cyanobacteria')
 add_numbered('Figure2','Taxonomic and functional population structure','NMDS_ABC_og_with_Cyanobacteria')
 add_numbered('Figure3','Population-centroid distances','Centroid_ABC_og_with_Cyanobacteria')
 add_numbered('Figure4','Dominant functions of the symbiotic partners','Partner_functions_2x2_subRole')
 add_numbered('Figure5',paste('Exploratory differential abundance:',selected,'versus Styria'),'Volcano_extreme_with_Cyanobacteria')
 suppkeys<-c('Taxonomy_Family_without_Cyanobacteria','NMDS_ABC_og_without_Cyanobacteria','NMDS_ABC_subRole_with_Cyanobacteria','NMDS_ABC_subRole_without_Cyanobacteria','Centroid_ABC_og_without_Cyanobacteria','Centroid_ABC_subRole_with_Cyanobacteria','Centroid_ABC_subRole_without_Cyanobacteria','Volcano_extreme_without_Cyanobacteria','Core_family_recovery','Core_ASV_detection','Alpha_family_boxplots_with_without_Cyanobacteria','Partner_functions_2x2_og','Shared_family_protein_contributions','Shared_family_functional_structure','Mantel_test_summary','Mantel_scatter_All_bacteria_all_seven','Mantel_scatter_Noncyano_bacteria_all_seven','Mantel_scatter_All_bacteria_Kilimanjaro','Mantel_scatter_Noncyano_bacteria_Kilimanjaro','QC_sequencing_depth','Core_definition_sensitivity')
 supptitles<-c('Family composition excluding Cyanobacteria','OG NMDS excluding Cyanobacteria','SubRole NMDS including Cyanobacteria','SubRole NMDS excluding Cyanobacteria','OG centroid distances excluding Cyanobacteria','SubRole centroid distances including Cyanobacteria','SubRole centroid distances excluding Cyanobacteria',paste('Exploratory',selected,'versus Styria contrasts excluding Cyanobacteria'),'Strict-core family recovery','Detection of abundant strict-core ASVs','Family-level alpha diversity with and without Cyanobacteria','Dominant annotated partner OGs','Shared-family protein contributions','Within-family functional differentiation','Mantel test summary','All-bacterial Mantel distances: seven populations','Non-cyanobacterial Mantel distances: seven populations','All-bacterial Mantel distances: Kilimanjaro populations','Non-cyanobacterial Mantel distances: Kilimanjaro populations','Sequencing depth','Core-definition and rarefaction sensitivity')
 stopifnot(length(suppkeys)==length(supptitles),!anyDuplicated(suppkeys))
 for(i in seq_along(suppkeys))add_numbered(paste0('Supp-F',i),supptitles[i],suppkeys[i])
 stopifnot(identical(names(numbered)[1:5],paste0('Figure',1:5)),length(numbered)==26L)
 manifest<-list();captions<-list()
 for(id in names(numbered)){
  entry<-numbered[[id]];section<-if(startsWith(id,'Figure'))'Main' else 'Supplementary';filename<-paste0(id,'.jpg')
  message('Exporting ',id,': ',entry$title)
  ggsave(file.path(fdir,section,filename),plot=entry$plot,device='jpeg',width=entry$width,height=entry$height,units='in',dpi=600,quality=100,bg='white',limitsize=TRUE)
  manifest[[id]]<-data.frame(figure=id,title=entry$title,relative_file=file.path(section,filename),source_plots=entry$source,width_in=entry$width,height_in=entry$height,dpi=600)
  captions[[id]]<-c(paste0(id,'. ',entry$title),entry$caption,'')
 }
 writeLines(unlist(captions),file.path(fdir,'FIGURE_CAPTIONS.txt'))
 write.csv(do.call(rbind,manifest),file.path(fdir,'FIGURE_INDEX.csv'),row.names=FALSE)
 writeLines(c('MAIN FIGURES: Figure1.jpg to Figure5.jpg in Main.','SUPPLEMENTARY FIGURES: Supp-F1.jpg to Supp-F21.jpg in Supplementary.','All exports are JPG, 600 dpi. Other figure directories are not modified.','Captions and figure titles: FIGURE_CAPTIONS.txt. Output mapping: FIGURE_INDEX.csv.','No class/order taxonomy or exploratory family-correlation figures are exported.','Completed primary native results are reused; the added whole-metaproteome exclusion model is fitted for panel C of Supp-F8.','Please inspect all rendered figures for clipping before manuscript submission.'),file.path(fdir,'READ_ME.txt'))
 cat('\nNumbered JPG set complete: 5 main figures and 21 supplementary figures.\nOutput: ',fdir,'\n',sep='')

}, envir=new.env(parent=context))

run_full <- function() {
 args <- commandArgs(trailingOnly=TRUE)
 root <- if(length(args)) normalizePath(args[1],mustWork=TRUE) else normalizePath(getwd(),mustWork=TRUE)
 if(!dir.exists(file.path(root,'input'))) stop('Set your working directory to Lobaria_reanalysis, or pass that directory as the first argument.')
 needed <- c('vegan','limma','statmod','DESeq2','ggplot2','ggrepel')
 missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
 if(length(missing)) stop('Missing packages: ',paste(missing,collapse=', '),'. Install CRAN packages vegan, ggplot2, ggrepel, statmod; install limma and DESeq2 using BiocManager.')
 out <- file.path(root,'R_results');dir.create(out,showWarnings=FALSE);dir.create(file.path(out,'matrices'),showWarnings=FALSE);dir.create(file.path(out,'figures'),showWarnings=FALSE);dir.create(file.path(out,'DA'),showWarnings=FALSE)
 logcon <- file(file.path(out,'console_log.txt'),open='wt');sink(logcon,split=TRUE);sink(logcon,type='message')
 on.exit({sink(type='message');sink();close(logcon)},add=TRUE)
 options(warn=1,width=160);set.seed(20261031)
 cat('LOBARIA R REANALYSIS\nStarted: ',format(Sys.time()),'\nProject: ',root,'\n',sep='')
 cat('No old result table is used to build abundance profiles. Shared draw indices are used when supplied.\n')
 rd <- function(n,...) read.csv(file.path(root,'input',n),check.names=FALSE,stringsAsFactors=FALSE,...)
 meta<-rd('meta_data_TZ_AU.csv');pops<-c('Mweka1','Mweka2','Mweka3','Marangu','Maua','Umbwe','Styria');meta<-meta[order(match(meta$group,pops),meta$replicate),];samples<-meta$id;g<-factor(meta$group,levels=pops);tz<-which(g!='Styria');au<-which(g=='Styria')
 C<-as.matrix(rd('amplicon_16S_asv_TZ_AU.csv',row.names=1)[,samples]);storage.mode(C)<-'double';tax<-rd('amplicon_16S_tax_TZ_AU.csv',row.names=1);tax<-tax[rownames(C),,drop=FALSE];tax[is.na(tax)]<-'unclassified';prot<-rd('metaproteome_TZ_AU.csv');prot[is.na(prot)]<-'unclassified';V<-as.matrix(prot[,samples]);storage.mode(V)<-'double'
 stopifnot(length(samples)==21L,!anyDuplicated(samples),all(table(g)==3),all(C>=0),all(C==round(C)),all(is.finite(V)),all(V>=0),all(abs(colSums(V)-1)<1e-5))
 org<-apply(tax,1,function(z)any(grepl('chloroplast|mitochond',z,ignore.case=TRUE)));keep<-tax$Kingdom=='Bacteria' & !org;C<-C[keep,,drop=FALSE];tax<-tax[keep,,drop=FALSE]
 prot$partner<-'Other';prot$partner[prot$superKingdom=='Bacteria' & prot$Phylum!='Cyanobacteria']<-'Noncyano_bacteria';prot$partner[prot$Phylum=='Ascomycota']<-'Mycobiont';prot$partner[prot$Phylum=='Chlorophyta']<-'Phycobiont';prot$partner[prot$superKingdom=='Bacteria' & prot$Phylum=='Cyanobacteria']<-'Cyanobiont'
 unknown<-c('unclassified','various','Function unknown','poorly characterized','unknown','NA','');named<-function(x)!x %in% unknown
 finiteBH<-function(p){z<-rep(NA_real_,length(p));ok<-is.finite(p);z[ok]<-p.adjust(p[ok],'BH');z}
 write<-function(x,n) {write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='');invisible(x)}
 relative<-function(x){stopifnot(all(rowSums(x)>0));sweep(x,1,rowSums(x),'/')}
 agg<-function(mask,le)t(rowsum(V[mask,,drop=FALSE],group=prot[[le]][mask],reorder=TRUE))
 core<-rowSums(C[,tz,drop=FALSE]>0)==18;shared<-core & rowSums(C[,au,drop=FALSE]>0)==3;families<-sort(unique(tax$Family[shared & named(tax$Family)]));ca<-data.frame(ASV=rownames(C)[core],tax[core,,drop=FALSE],Styria_detections=rowSums(C[core,au,drop=FALSE]>0),mean_Kilimanjaro_fraction=rowMeans(sweep(C[core,tz,drop=FALSE],2,colSums(C[,tz,drop=FALSE]),'/')),mean_Styria_fraction=rowMeans(sweep(C[core,au,drop=FALSE],2,colSums(C[,au,drop=FALSE]),'/')))
 write(ca,'core_ASVs');threshold<-list()
 for(co in c('All_bacteria','Noncyano_bacteria')) {
  A<-C[if(co=='All_bacteria')rep(TRUE,nrow(C)) else tax$Phylum!='Cyanobacteria',,drop=FALSE];rr<-sweep(A,2,colSums(A),'/')
  for(def in c('one_read','two_reads','relative_0.01pct')) for(n in c(18,16,14)) {
   det<-switch(def,one_read=A>0,two_reads=A>=2,relative_0.01pct=rr>=.0001);k<-rowSums(det[,tz,drop=FALSE])>=n
   threshold[[length(threshold)+1]]<-data.frame(compartment=co,definition=def,min_Kilimanjaro_detections=n,core_ASVs=sum(k),Styria_any=sum(rowSums(det[k,au,drop=FALSE])>0),Styria_all3=sum(rowSums(det[k,au,drop=FALSE])==3))
  }
 }
 write(do.call(rbind,threshold),'core_thresholds');depth<-min(colSums(C));rare<-list()
 for(b in 1:200) {Z<-t(vegan::rrarefy(t(C),sample=depth));k<-rowSums(Z[,tz,drop=FALSE]>0)==18;rare[[b]]<-data.frame(draw=b,depth=depth,core_ASVs=sum(k),Styria_any=sum(rowSums(Z[k,au,drop=FALSE]>0)>0),Styria_all3=sum(rowSums(Z[k,au,drop=FALSE]>0)==3))}
 write(do.call(rbind,rare),'core_rarefaction');write(data.frame(clean_ASVs=nrow(C),core_ASVs=sum(core),Styria_any=sum(ca$Styria_detections>0),Styria_all3=sum(shared),mean_core_Kilimanjaro=sum(ca$mean_Kilimanjaro_fraction),mean_core_Styria=sum(ca$mean_Styria_fraction),shared_named_families=length(families),rarefaction_depth=depth),'core_summary');cat('\nCORE SUMMARY\n');print(read.csv(file.path(out,'core_summary.csv')))
 write(data.frame(meta,clean_bacterial_reads=colSums(C),detected_protein_groups=colSums(V>0),total_NSAF=colSums(V)),'sample_QC')
 masks<-list(Whole=rep(TRUE,nrow(prot)),All_bacteria=prot$superKingdom=='Bacteria',Noncyano_bacteria=prot$partner=='Noncyano_bacteria',Cyanobiont=prot$partner=='Cyanobiont',Mycobiont=prot$partner=='Mycobiont',Phycobiont=prot$partner=='Phycobiont',Shared_core_lineages=prot$partner=='Noncyano_bacteria' & prot$Family %in% families)
 profiles<-list();counts<-list();keys<-list()
 add<-function(data,co,le,x){tag<-paste(data,co,le,sep='__');profiles[[tag]]<<-relative(x);keys[[tag]]<<-c(data,co,le)}
 for(co in c('All_bacteria','Noncyano_bacteria')) for(le in c('Class','Order','Family')) {
  sel<-if(co=='All_bacteria')rep(TRUE,nrow(C)) else tax$Phylum!='Cyanobacteria';x<-t(rowsum(C[sel,,drop=FALSE],tax[[le]][sel],reorder=TRUE));x<-x[,colSums(x)>=20 & colSums(x>0)>=3,drop=FALSE];counts[[paste(co,le,sep='__')]]<-x;add('16S',co,le,x);add('Protein_taxonomy',co,le,agg(masks[[co]],le))
 }
 for(co in names(masks))for(le in c('og','subRole','mainRole'))add('Function',co,le,agg(masks[[co]],le))
 B<-4999;NB<-1999
 drawfile<-function(n)file.path(root,'python_results',n)
 perms<-if(file.exists(drawfile('permutations_1based.csv')))as.matrix(read.csv(drawfile('permutations_1based.csv'),header=FALSE)) else t(replicate(B,sample.int(21)))
 boot<-if(file.exists(drawfile('bootstrap_1based.csv')))as.matrix(read.csv(drawfile('bootstrap_1based.csv'),header=FALSE)) else t(replicate(NB,unlist(lapply(pops,function(z)sample(which(g==z),3,replace=TRUE)))))
 stopifnot(ncol(perms)==21,ncol(boot)==21);B<-nrow(perms);NB<-nrow(boot)
 # Positive full-space PCoA coordinates with Lingoes correction; same definition as Python.
 pcoa<-function(D){n<-nrow(D);H<-diag(n)-matrix(1/n,n,n);E<-eigen(-.5*H%*%(D^2)%*%H,symmetric=TRUE);cc<-max(0,-min(E$values));Ec<-eigen(-.5*H%*%(D^2+2*cc*(1-diag(n)))%*%H,symmetric=TRUE);k<-Ec$values>1e-10;list(Z=sweep(Ec$vectors[,k,drop=FALSE],2,sqrt(Ec$values[k]),'*'),constant=cc)}
 pm<-function(D,gg,pr){dd<-as.dist(D);ans<-vegan::adonis2(dd~gg,permutations=pr);c(F=ans$F[1],R2=ans$R2[1],p=ans$`Pr(>F)`[1])}
 # Label-permutation dispersion test: recompute group centroids for every permutation.
 dispersion<-function(Z,gg,pr){lev<-unique(gg);n<-length(gg);stat<-function(ids){zz<-Z[ids,,drop=FALSE];ds<-numeric(n);for(a in lev){ii<-which(gg==a);ds[ii]<-sqrt(rowSums(sweep(zz[ii,,drop=FALSE],2,colMeans(zz[ii,,drop=FALSE]),'-')^2))*sqrt(length(ii)/(length(ii)-1))};means<-sapply(lev,function(a)mean(ds[gg==a]));ssb<-sum(sapply(seq_along(lev),function(j)sum(gg==lev[j])*(means[j]-mean(ds))^2));ssw<-sum(sapply(seq_along(lev),function(j)sum((ds[gg==lev[j]]-means[j])^2)));(ssb/(length(lev)-1))/(ssw/(n-length(lev)))};f<-stat(seq_len(n));c(F=f,p=(1+sum(apply(pr,1,stat)>=f-1e-12))/(nrow(pr)+1))}
 cache<-list();cent<-list();rboots<-list();global<-list();centrows<-list();sites<-list()
 for(tag in names(profiles)) {
  x<-profiles[[tag]];k<-keys[[tag]];D<-as.matrix(vegan::vegdist(x,method='bray'));pc<-pcoa(D);Z<-pc$Z;CT<-t(vapply(pops,function(z)colMeans(Z[g==z,,drop=FALSE]),numeric(ncol(Z))));CD<-as.matrix(dist(CT));dimnames(CD)<-list(pops,pops);cache[[tag]]<-D;cent[[tag]]<-CD
  write.csv(data.frame(sample=samples,x,check.names=FALSE),file.path(out,'matrices',paste0(tag,'.csv')),row.names=FALSE);write.csv(CD,file.path(out,'matrices',paste0('centroids__',tag,'.csv')))
  q<-pm(D,g,perms);di<-dispersion(Z,g,perms);global[[length(global)+1]]<-data.frame(dataset=k[1],compartment=k[2],level=k[3],F=q['F'],R2=q['R2'],p=q['p'],dispersion_F=di['F'],dispersion_p=di['p'])
  MP<-as.matrix(vegan::vegdist(t(vapply(pops,function(z)colMeans(x[g==z,,drop=FALSE]),numeric(ncol(x)))),method='bray'))
  methods<-list(full_PCoA_Lingoes=CD,meanprofile_Bray_sensitivity=MP)
  for(method in names(methods)) {
   DD<-methods[[method]];within<-DD[1:6,1:6][upper.tri(DD[1:6,1:6])];st<-DD[1:6,7];difference<-mean(st)-mean(within)
   refs<-sapply(1:7,function(r){oth<-setdiff(1:7,r);ww<-DD[oth,oth];mean(DD[r,oth])-mean(ww[upper.tri(ww)])})
   row<-data.frame(dataset=k[1],compartment=k[2],level=k[3],method=method,mean_within=mean(within),mean_Styria=mean(st),difference=difference,ratio=mean(st)/mean(within),reference_exact_p=mean(refs>=difference-1e-12),all_five_count=sum(sapply(1:6,function(i)st[i]>max(DD[i,setdiff(1:6,i)]))),lingoes_constant=pc$constant,bootstrap_low=NA_real_,bootstrap_high=NA_real_)
   if(method=='full_PCoA_Lingoes') {bv<-apply(boot,1,function(ii){zc<-t(vapply(pops,function(z)colMeans(Z[ii[g==z],,drop=FALSE]),numeric(ncol(Z))));dd<-as.matrix(dist(zc));ww<-dd[1:6,1:6];mean(dd[1:6,7])-mean(ww[upper.tri(ww)])});row$bootstrap_low<-quantile(bv,.025);row$bootstrap_high<-quantile(bv,.975)}
   centrows[[length(centrows)+1]]<-row
   for(i in 1:6){ds<-DD[i,setdiff(1:6,i)];sites[[length(sites)+1]]<-data.frame(dataset=k[1],compartment=k[2],level=k[3],method=method,population=pops[i],distance_Styria=st[i],mean_other_Kilimanjaro=mean(ds),min_other_Kilimanjaro=min(ds),max_other_Kilimanjaro=max(ds),n_other_Kilimanjaro_closer=sum(ds<st[i]),Styria_farther_than_all_five=st[i]>max(ds))}
  }
  rboots[[tag]]<-apply(boot,1,function(ii){dd<-D[ii,ii]^2;tot<-sum(dd)/42;wi<-sum(sapply(pops,function(z)sum(dd[g==z,g==z])/6));(tot-wi)/tot})
  cat('Completed population profile:',tag,' R2=',q['R2'],' P=',q['p'],'\n')
 }
 global<-do.call(rbind,global);global$q_population<-p.adjust(global$p,'BH');global$q_dispersion<-p.adjust(global$dispersion_p,'BH');write(global,'population_PERMANOVA');centrows<-do.call(rbind,centrows);write(centrows,'centroid_summary');sites<-do.call(rbind,sites);write(sites,'centroid_each_population')
 pr18<-if(file.exists(drawfile('permutations18_1based.csv')))as.matrix(read.csv(drawfile('permutations18_1based.csv'),header=FALSE)) else t(replicate(B,sample.int(18)))
 wk<-list();resolution<-list()
 for(co in c('Whole','All_bacteria','Noncyano_bacteria','Cyanobiont','Mycobiont','Phycobiont')) {
  for(le in c('og','subRole','mainRole')){tag<-paste('Function',co,le,sep='__');z<-pm(cache[[tag]][tz,tz],droplevels(g[tz]),pr18);wk[[length(wk)+1]]<-data.frame(compartment=co,level=le,F=z['F'],R2=z['R2'],p=z['p'])}
  for(pair in list(c('og','subRole'),c('subRole','mainRole'),c('og','mainRole'))){v<-rboots[[paste('Function',co,pair[2],sep='__')]]-rboots[[paste('Function',co,pair[1],sep='__')]];resolution[[length(resolution)+1]]<-data.frame(compartment=co,from_level=pair[1],to_level=pair[2],bootstrap_low=quantile(v,.025),bootstrap_high=quantile(v,.975))}
 }
 wk<-do.call(rbind,wk);wk$q<-p.adjust(wk$p,'BH');write(wk,'within_Kilimanjaro_PERMANOVA');write(do.call(rbind,resolution),'resolution_paired_bootstrap')
 # Explicit two-sided Mantel test (vegan::mantel's usual alternative is positive).
 mantel2<-function(A,D,pr){ij<-which(upper.tri(A),arr.ind=TRUE);a<-A[ij];r<-cor(a,D[ij]);rp<-apply(pr,1,function(ii)cor(a,D[ii,ii][ij]));c(r=r,p=(1+sum(abs(rp)>=abs(r)-1e-12))/(nrow(pr)+1))}
 mm<-list()
 for(scope in c('all_seven','Kilimanjaro')){ix<-if(scope=='all_seven')1:21 else tz;pr<-if(scope=='all_seven')perms else pr18
  for(co in c('All_bacteria','Noncyano_bacteria'))for(tl in c('Class','Order','Family'))for(fl in c('og','subRole','mainRole')){z<-mantel2(cache[[paste('16S',co,tl,sep='__')]][ix,ix],cache[[paste('Function',co,fl,sep='__')]][ix,ix],pr);mm[[length(mm)+1]]<-data.frame(scope=scope,compartment=co,taxonomy_level=tl,function_level=fl,r=z['r'],p=z['p'])}}
 mm<-do.call(rbind,mm);mm$q<-ave(mm$p,mm$scope,FUN=function(p)p.adjust(p,'BH'));write(mm,'Mantel')
 # Dominant functions, within-partner denominator always includes unknown annotations.
 rankings<-list();by_pop<-list();shares<-list()
 for(co in c('Cyanobiont','Mycobiont','Noncyano_bacteria','Phycobiont')) {
  mass<-colSums(V[masks[[co]],,drop=FALSE]);shares[[co]]<-data.frame(sample=samples,population=as.character(g),partner=co,whole_metaproteome_share=mass)
  for(le in c('og','subRole')){x<-profiles[[paste('Function',co,le,sep='__')]];mu<-colMeans(x);rr<-rank(-mu,ties.method='min');nr<-rep(NA_real_,length(mu));nr[named(names(mu))]<-rank(-mu[named(names(mu))],ties.method='min')
   desc<-if(le=='og')vapply(colnames(x),function(z)paste(unique(prot$desc[masks[[co]] & prot[[le]]==z]),collapse=' | '),character(1)) else colnames(x)
   rankings[[length(rankings)+1]]<-data.frame(partner=co,level=le,feature=colnames(x),description=desc,mean_fraction=mu,sd_fraction=apply(x,2,sd),rank_all=rr,rank_named=nr,n_thalli_detected=colSums(x>0),n_populations_detected_2of3=rowSums(vapply(pops,function(z)colSums(x[g==z,,drop=FALSE]>0)>=2,logical(ncol(x)))),named=named(colnames(x)))
   for(pop in pops){xx<-x[g==pop,,drop=FALSE];by_pop[[length(by_pop)+1]]<-data.frame(partner=co,level=le,population=pop,feature=colnames(x),mean_within_partner=colMeans(xx),sd_within_partner=apply(xx,2,sd),n_detected=colSums(xx>0),mean_whole_metaproteome=colMeans(sweep(xx,1,mass[g==pop],'*')))}
  }
 }
 rankings<-do.call(rbind,rankings);rankings$informative_annotation<-named(rankings$description) & !grepl('unknown function|domain of unknown|DUF[0-9]',rankings$description,ignore.case=TRUE);rankings$rank_informative<-NA_real_;ii<-which(rankings$informative_annotation);rankings$rank_informative[ii]<-ave(-rankings$mean_fraction[ii],interaction(rankings$partner[ii],rankings$level[ii]),FUN=function(x)rank(x,ties.method='min'));write(rankings,'partner_function_rankings');write(do.call(rbind,by_pop),'partner_functions_by_population');write(do.call(rbind,shares),'partner_whole_metaproteome_shares')
 # Family-specific functional profiles and representation are distinct quantities.
 coverage<-list();famtests<-list();fammean<-list();family_coeff<-list();family_for_limma<-matrix(0,length(families),21,dimnames=list(families,samples))
 for(fam in families) {
  mask<-prot$partner=='Noncyano_bacteria' & prot$Family==fam;mass<-colSums(V[mask,,drop=FALSE]);nd<-sapply(pops,function(z)sum(mass[g==z]>0));bden<-colSums(V[masks$All_bacteria,,drop=FALSE]);nden<-colSums(V[masks$Noncyano_bacteria,,drop=FALSE]);family_for_limma[fam,]<-mass/bden
  coverage[[fam]]<-data.frame(family=fam,sample=samples,population=as.character(g),protein_groups=sum(mask),detected_protein_groups=colSums(V[mask,,drop=FALSE]>0),whole_fraction=mass,noncyano_fraction=mass/nden,all_bacterial_fraction=mass/bden)
  for(le in c('og','subRole')) {
   if(!any(mask))next
   x<-agg(mask,le);active<-rowSums(x)>0;y<-relative(x[active,,drop=FALSE]);gg<-droplevels(g[active]);write.csv(data.frame(sample=samples[active],y,check.names=FALSE),file.path(out,'matrices',paste0('family__',fam,'__',le,'.csv')),row.names=FALSE)
   for(pop in pops){sel<-gg==pop;if(!any(sel))next;fammean[[length(fammean)+1]]<-data.frame(family=fam,level=le,population=pop,feature=colnames(y),mean_within_family=colMeans(y[sel,,drop=FALSE]),n_positive_thalli=sum(sel))}
   eligible<-all(nd>=2) && ncol(y)>=2 && any(apply(y,2,function(v)diff(range(v)))>0)
   row<-data.frame(family=fam,level=le,protein_groups=sum(mask),n_positive_thalli=nrow(y),min_positive_per_population=min(nd),eligible_all_seven=eligible,F=NA_real_,R2=NA_real_,p=NA_real_,dispersion_p=NA_real_,mean_Styria_distance=NA_real_,mean_Kilimanjaro_distance=NA_real_)
   if(eligible){D<-as.matrix(vegan::vegdist(y));pf<-drawfile(paste0('family_permutations__',fam,'__',le,'.csv'));pr<-if(file.exists(pf))as.matrix(read.csv(pf,header=FALSE)) else if(all(active))perms else t(replicate(B,sample.int(nrow(y))));z<-pm(D,gg,pr);pc<-pcoa(D);di<-dispersion(pc$Z,gg,pr);ct<-t(vapply(pops,function(pop)colMeans(pc$Z[gg==pop,,drop=FALSE]),numeric(ncol(pc$Z))));cd<-as.matrix(dist(ct));ww<-cd[1:6,1:6];row$F<-z['F'];row$R2<-z['R2'];row$p<-z['p'];row$dispersion_p<-di['p'];row$mean_Styria_distance<-mean(cd[1:6,7]);row$mean_Kilimanjaro_distance<-mean(ww[upper.tri(ww)])}
   famtests[[length(famtests)+1]]<-row
  }
  for(pair in combn(pops,2,simplify=FALSE)){a<-family_for_limma[fam,g==pair[1]];b<-family_for_limma[fam,g==pair[2]];if(sum(a>0)>=2 && sum(b>0)>=2)family_coeff[[length(family_coeff)+1]]<-data.frame(family=fam,popA=pair[1],popB=pair[2],log2FoldChange=mean(log2(b[b>0]))-mean(log2(a[a>0])),nA=sum(a>0),nB=sum(b>0))}
 }
 write(do.call(rbind,coverage),'shared_family_contributions');write(do.call(rbind,fammean),'shared_family_functions_by_population');ft<-do.call(rbind,famtests);ft$q<-finiteBH(ft$p);write(ft,'shared_family_profile_tests');write(do.call(rbind,family_coeff),'shared_family_detected_log_coefficients')
 # Native R moderation for family contribution; zeros are nondetections, not imputed.
 Y<-log2(family_for_limma);Y[!is.finite(Y)]<-NA;design<-model.matrix(~0+g);colnames(design)<-pops;nd<-sapply(pops,function(z)rowSums(is.finite(Y[,g==z,drop=FALSE])));base<-rowSums(is.finite(Y))>=6 & rowSums(is.finite(Y))-rowSums(nd>0)>=3
 if(sum(base)>=2){fit<-limma::lmFit(Y[base,,drop=FALSE],design);ans<-list()
  for(pair in combn(pops,2,simplify=FALSE)){ids<-rownames(Y)[base & nd[,pair[1]]>=2 & nd[,pair[2]]>=2];if(length(ids)<2)next;ct<-matrix(0,7,1,dimnames=list(pops,NULL));ct[pair[1],1]<- -1;ct[pair[2],1]<-1;ff<-limma::eBayes(limma::contrasts.fit(fit[match(ids,rownames(fit)),],ct),trend=FALSE,robust=FALSE);tt<-limma::topTable(ff,number=Inf,sort.by='none');ans[[length(ans)+1]]<-data.frame(family=rownames(tt),popA=pair[1],popB=pair[2],log2FoldChange=tt$logFC,p=tt$P.Value)}
  if(length(ans)){aa<-do.call(rbind,ans);aa$q_across_families_and_21_pairs<-p.adjust(aa$p,'BH');write(aa,'shared_family_contribution_limma')}
 }
 # Graphics exported after all native feature models finish.
 
 # Embedded module: secondary_R.R
 eval(quote({
# Secondary analyses; source within main() after primary matrices are available.
ss<-list()
for(co in c('Whole','All_bacteria','Noncyano_bacteria','Cyanobiont','Mycobiont','Phycobiont'))for(le in c('og','subRole','mainRole'))for(mode in c('exclude_mw23','remove_uninformative_bins')){
 x<-profiles[[paste('Function',co,le,sep='__')]];active<-if(mode=='exclude_mw23')samples!='mw23' else rep(TRUE,length(samples));if(mode=='remove_uninformative_bins')x<-x[,named(colnames(x)),drop=FALSE];y<-relative(x[active,,drop=FALSE]);gg<-g[active];Z<-pcoa(as.matrix(vegan::vegdist(y)))$Z;CT<-t(vapply(pops,function(pop)colMeans(Z[gg==pop,,drop=FALSE]),numeric(ncol(Z))));D<-as.matrix(dist(CT));ww<-D[1:6,1:6];st<-D[1:6,7];ss[[length(ss)+1]]<-data.frame(compartment=co,level=le,sensitivity=mode,mean_within=mean(ww[upper.tri(ww)]),mean_Styria=mean(st),ratio=mean(st)/mean(ww[upper.tri(ww)]),all_five_count=sum(sapply(1:6,function(i)st[i]>max(D[i,setdiff(1:6,i)]))))
}
write(do.call(rbind,ss),'centroid_sensitivities')
allperms<-function(v){if(length(v)==1)return(matrix(v,1,1));do.call(rbind,lapply(seq_along(v),function(i)cbind(v[i],allperms(v[-i]))))}
ms<-list()
for(scope in c('all_seven','Kilimanjaro')){
 n<-if(scope=='all_seven')21 else 18;nc<-n/3;ix<-seq_len(n);pf<-drawfile(paste0('Mantel_within_population_perms_',n,'_1based.csv'));pr<-if(file.exists(pf))as.matrix(read.csv(pf,header=FALSE)) else t(replicate(B,unlist(lapply(seq_len(nc),function(j)sample((3*j-2):(3*j))))));pc<-allperms(seq_len(nc))
 for(co in c('All_bacteria','Noncyano_bacteria'))for(tl in c('Class','Order','Family'))for(fl in c('og','subRole','mainRole'))for(method in c('within_population','population_centroids')){
  if(method=='within_population'){A<-cache[[paste('16S',co,tl,sep='__')]][ix,ix];D<-cache[[paste('Function',co,fl,sep='__')]][ix,ix];z<-mantel2(A,D,pr)}else{A<-cent[[paste('16S',co,tl,sep='__')]][seq_len(nc),seq_len(nc)];D<-cent[[paste('Function',co,fl,sep='__')]][seq_len(nc),seq_len(nc)];ij<-which(upper.tri(A),arr.ind=TRUE);r<-cor(A[ij],D[ij]);rp<-apply(pc,1,function(ii)cor(A[ij],D[ii,ii][ij]));z<-c(r=r,p=mean(abs(rp)>=abs(r)-1e-12))}
  ms[[length(ms)+1]]<-data.frame(scope=scope,compartment=co,taxonomy_level=tl,function_level=fl,method=method,r=z['r'],p=z['p'])
 }
}
ms<-do.call(rbind,ms);ms$q<-ave(ms$p,interaction(ms$scope,ms$method),FUN=finiteBH);write(ms,'Mantel_secondary')

 }), envir=environment())

 # Cross-language checks of deterministic quantities, with explicit tolerances.
 compare<-function(name,keys,values,tolerance=1e-8){py<-file.path(root,'python_results',paste0(name,'.csv'));rr<-file.path(out,paste0(name,'.csv'));if(!file.exists(py)||!file.exists(rr))return(NULL);a<-read.csv(py,check.names=FALSE);b<-read.csv(rr,check.names=FALSE);if(name=='shared_family_profile_tests' && 'record_type' %in% names(a))a<-a[a$record_type=='global_test',];m<-merge(a,b,by=keys,suffixes=c('_python','_R'));do.call(rbind,lapply(values,function(v){aa<-m[[paste0(v,'_python')]];bb<-m[[paste0(v,'_R')]];ok<-is.finite(aa)&is.finite(bb);err<-if(any(ok))max(abs(aa[ok]-bb[ok])) else NA_real_;data.frame(table=name,quantity=v,n_compared=sum(ok),max_absolute_error=err,tolerance=tolerance,passed=if(any(ok))err<=tolerance else NA)}))}
 checks<-rbind(compare('centroid_summary',c('dataset','compartment','level','method'),c('mean_within','mean_Styria','difference','ratio','reference_exact_p')),compare('population_PERMANOVA',c('dataset','compartment','level'),c('F','R2','p')),compare('Mantel',c('scope','compartment','taxonomy_level','function_level'),c('r','p')),compare('partner_function_rankings',c('partner','level','feature'),c('mean_fraction','sd_fraction')),compare('shared_family_detected_log_coefficients',c('family','popA','popB'),c('log2FoldChange')),compare('centroid_sensitivities',c('compartment','level','sensitivity'),c('mean_within','mean_Styria','ratio')),compare('Mantel_secondary',c('scope','compartment','taxonomy_level','function_level','method'),c('r','p')),compare('shared_family_profile_tests',c('family','level'),c('F','R2','p')))
 if(!is.null(checks)){write(checks,'Python_R_checks');cat('\nCROSS-LANGUAGE CHECKS\n');print(checks)}
 # Run full original native limma analyses and native DESeq2 models.
 cat('\nRUNNING NATIVE R FEATURE MODELS\n')
 
 # Embedded module: limma_functions.R
 eval(quote({
# R 4.3.3, limma 3.58.1, statmod 1.5.0.
# Seven-population models of log2 detected relative NSAF; no imputation.
suppressPackageStartupMessages(library(limma))
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1]) else normalizePath('.')
out <- file.path(root,'R_results','DA')
meta <- read.csv(file.path(root,'input','meta_data_TZ_AU.csv'),check.names=FALSE)
prot <- read.csv(file.path(root,'input','metaproteome_TZ_AU.csv'),check.names=FALSE)
tax <- read.csv(file.path(root,'input','amplicon_16S_tax_TZ_AU.csv'),row.names=1,check.names=FALSE)
asv <- read.csv(file.path(root,'input','amplicon_16S_asv_TZ_AU.csv'),row.names=1,check.names=FALSE)
samples <- meta$id;pops <- unique(meta$group);g <- factor(meta$group,levels=pops);design <- model.matrix(~0+g);colnames(design)<-pops
prot$partner <- 'Other'
prot$partner[prot$superKingdom=='Bacteria' & prot$Phylum!='Cyanobacteria'] <- 'Noncyano_bacteria'
prot$partner[prot$Phylum=='Ascomycota'] <- 'Mycobiont'
prot$partner[prot$Phylum=='Chlorophyta'] <- 'Phycobiont'
prot$partner[prot$superKingdom=='Bacteria' & prot$Phylum=='Cyanobacteria'] <- 'Cyanobiont'
org <- apply(tax,1,function(z) any(grepl('chloroplast|mitochond',z,ignore.case=TRUE)))
core <- rowSums(asv[,samples]>0)==21 & tax$Kingdom=='Bacteria' & !org
families <- unique(na.omit(tax$Family[core]));families<-setdiff(families,c('unclassified','various'))
masks <- list(Whole=rep(TRUE,nrow(prot)),All_bacteria=prot$superKingdom=='Bacteria',Noncyano_bacteria=prot$partner=='Noncyano_bacteria',Cyanobiont=prot$partner=='Cyanobiont',Phycobiont=prot$partner=='Phycobiont',Mycobiont=prot$partner=='Mycobiont',Shared_core_lineages=prot$partner=='Noncyano_bacteria' & prot$Family %in% families)
stopifnot(sum(masks$Shared_core_lineages)==259)
pair_results<-list();overall<-list();det_results<-list();diag_results<-list();cc_results<-list();counter<-0
moderate <- function(fit) {
 n<-nrow(fit$coefficients)
 # Broad levels have too few features to justify a fitted intensity trend.
 trend <- n>=30; robust <- n>=20
 ans<-tryCatch(eBayes(fit,trend=trend,robust=robust),error=function(e){warning(paste('Empirical-Bayes fallback to no trend/nonrobust:',conditionMessage(e)));eBayes(fit,trend=FALSE,robust=FALSE)})
 ans
}
for(comp in names(masks)) for(lev in c('og','subRole','mainRole')) {
 mask<-masks[[comp]];v<-as.matrix(prot[mask,samples]);storage.mode(v)<-'double'
 M<-rowsum(v,group=prot[[lev]][mask],reorder=FALSE);M<-sweep(M,2,colSums(v),'/')
 write.csv(data.frame(feature=rownames(M),M,check.names=FALSE),file.path(out,paste0('functional_matrix_',comp,'_',lev,'.csv')),row.names=FALSE)
 named<-!rownames(M) %in% c('unclassified','various','Function unknown','poorly characterized')
 M<-M[named,,drop=FALSE];Y<-log2(M);Y[!is.finite(Y)]<-NA
 nd<-sapply(pops,function(z) rowSums(M[,g==z,drop=FALSE]>0));rownames(nd)<-rownames(M)
 df<-rowSums(is.finite(Y))-rowSums(nd>0)
 basekeep<-rowSums(is.finite(Y))>=6 & df>=3
 if(sum(basekeep)<2) next
 fit<-lmFit(Y[basekeep,,drop=FALSE],design);baseids<-rownames(fit$coefficients)
 # Correct omnibus test is equality of the seven means, not seven means vs zero.
 glob<-apply(nd[baseids,,drop=FALSE]>=2,1,all)
 if(sum(glob)>=2) {
  C<-sapply(2:7,function(j) {z<-rep(0,7);z[j]<-1;z[1]<- -1;z});rownames(C)<-pops
  fg<-moderate(contrasts.fit(fit[glob,],C));tt<-topTableF(fg,number=Inf,sort.by='none',adjust.method='BH')
  overall[[length(overall)+1]]<-data.frame(compartment=comp,level=lev,feature=rownames(tt),F=tt$F,p=tt$P.Value,q_within_level=tt$adj.P.Val,stringsAsFactors=FALSE)
 }
 for(ii in combn(7,2,simplify=FALSE)) {
  a<-pops[ii[1]];b<-pops[ii[2]];na<-nd[,a];nb<-nd[,b];eligible<-na>=2 & nb>=2 & basekeep
  det_results[[length(det_results)+1]]<-data.frame(compartment=comp,level=lev,popA=a,popB=b,feature=rownames(M),n_detected_A=na,n_detected_B=nb,arithmetic_mean_A=rowMeans(M[,g==a,drop=FALSE]),arithmetic_mean_B=rowMeans(M[,g==b,drop=FALSE]),quantitative_eligible=eligible,on_off_3vs0=(na==3 & nb==0)|(na==0 & nb==3),stringsAsFactors=FALSE)
  ids<-rownames(M)[eligible];if(length(ids)<2) next
  contrast<-matrix(0,7,1,dimnames=list(pops,paste0(b,'_minus_',a)));contrast[a,1]<- -1;contrast[b,1]<-1
  fc<-moderate(contrasts.fit(fit[match(ids,baseids),],contrast));tt<-topTable(fc,number=Inf,sort.by='none',adjust.method='BH',confint=TRUE)
  pair_results[[length(pair_results)+1]]<-data.frame(compartment=comp,level=lev,popA=a,popB=b,feature=rownames(tt),log2FoldChange=tt$logFC,CI_L=tt$CI.L,CI_R=tt$CI.R,t=tt$t,pvalue=tt$P.Value,padj=tt$adj.P.Val,AveExpr=tt$AveExpr,n_detected_A=na[rownames(tt)],n_detected_B=nb[rownames(tt)],stringsAsFactors=FALSE)
  cc<-na[ids]==3 & nb[ids]==3
  if(sum(cc)>=2) {
   fcc<-moderate(contrasts.fit(fit[match(ids[cc],baseids),],contrast));tc<-topTable(fcc,number=Inf,sort.by='none',adjust.method='BH')
   cc_results[[length(cc_results)+1]]<-data.frame(compartment=comp,level=lev,popA=a,popB=b,feature=rownames(tc),log2FoldChange=tc$logFC,pvalue=tc$P.Value,padj=tc$adj.P.Val,stringsAsFactors=FALSE)
  }
 }
 diag_results[[length(diag_results)+1]]<-data.frame(compartment=comp,level=lev,n_named_features=nrow(M),n_base_model_features=length(baseids),n_global_all7_eligible=sum(glob),min_residual_df=min(fit$df.residual),max_residual_df=max(fit$df.residual))
 cat('Complete:',comp,lev,'model features',length(baseids),'\n');flush.console()
}
pair<-do.call(rbind,pair_results);pair$q_across_21_contrasts<-ave(pair$pvalue,interaction(pair$compartment,pair$level),FUN=function(p) p.adjust(p,'BH'))
write.csv(pair,file.path(out,'functional_pairwise_limma.csv'),row.names=FALSE)
cc<-do.call(rbind,cc_results);cc$q_across_21_contrasts<-ave(cc$pvalue,interaction(cc$compartment,cc$level),FUN=function(p) p.adjust(p,'BH'));write.csv(cc,file.path(out,'functional_pairwise_limma_completecase.csv'),row.names=FALSE)
ov<-do.call(rbind,overall);ov$q_across_functional_omnibus<-p.adjust(ov$p,'BH');write.csv(ov,file.path(out,'functional_omnibus_limma_F.csv'),row.names=FALSE)
write.csv(do.call(rbind,det_results),file.path(out,'functional_detection_results.csv'),row.names=FALSE)
write.csv(do.call(rbind,diag_results),file.path(out,'functional_model_diagnostics.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('All functional models complete\n')

 }), envir=new.env(parent=globalenv()))

 
 # Embedded module: DESeq2_taxonomy.R
 eval(quote({
# Native R DESeq2, same raw-count models/contrasts and BH family as Python PyDESeq2.
suppressPackageStartupMessages(library(DESeq2))
args<-commandArgs(trailingOnly=TRUE);root<-if(length(args))normalizePath(args[1]) else normalizePath(getwd())
out<-file.path(root,'R_results','DA');dir.create(out,recursive=TRUE,showWarnings=FALSE)
meta<-read.csv(file.path(root,'input','meta_data_TZ_AU.csv'),check.names=FALSE);samples<-meta$id;pops<-unique(meta$group)
tax<-read.csv(file.path(root,'input','amplicon_16S_tax_TZ_AU.csv'),row.names=1,check.names=FALSE);tax[is.na(tax)]<-'unclassified'
counts<-as.matrix(read.csv(file.path(root,'input','amplicon_16S_asv_TZ_AU.csv'),row.names=1,check.names=FALSE)[,samples]);tax<-tax[rownames(counts),,drop=FALSE]
org<-apply(tax,1,function(z)any(grepl('chloroplast|mitochond',z,ignore.case=TRUE)));bact<-tax$Kingdom=='Bacteria' & !org
md<-data.frame(group=factor(meta$group,levels=pops),row.names=samples);ans<-list();diag<-list()
for(co in c('All_bacteria','Noncyano_bacteria')) for(le in c('Family','Order','Class')) {
 mask<-bact & if(co=='All_bacteria')rep(TRUE,nrow(tax)) else tax$Phylum!='Cyanobacteria'
 M<-rowsum(counts[mask,,drop=FALSE],tax[[le]][mask],reorder=TRUE);M<-M[rowSums(M)>=20 & rowSums(M>0)>=3,,drop=FALSE];storage.mode(M)<-'integer'
 dds<-DESeqDataSetFromMatrix(M,colData=md,design=~group)
 dds<-DESeq(dds,fitType='parametric',sfType='poscounts',minReplicatesForReplace=7,quiet=TRUE)
 conv<-mcols(dds)$betaConv;if(is.null(conv))conv<-rep(TRUE,nrow(dds))
 for(pair in combn(pops,2,simplify=FALSE)) {
  a<-pair[1];b<-pair[2];rr<-as.data.frame(results(dds,contrast=c('group',b,a),independentFiltering=FALSE,cooksCutoff=TRUE,alpha=.05));rr$pvalue[!conv]<-NA_real_
  ans[[length(ans)+1]]<-data.frame(compartment=co,level=le,popA=a,popB=b,feature=rownames(rr),baseMean=rr$baseMean,log2FoldChange=rr$log2FoldChange,lfcSE=rr$lfcSE,stat=rr$stat,pvalue=rr$pvalue,n_detected_A=rowSums(M[,md$group==a,drop=FALSE]>0),n_detected_B=rowSums(M[,md$group==b,drop=FALSE]>0),LFC_converged=conv)
 }
 diag[[length(diag)+1]]<-data.frame(compartment=co,level=le,n_features=nrow(M),nonconverged=sum(!conv),fit_type=attr(dispersionFunction(dds),'fitType'));cat('Native DESeq2 complete:',co,le,'\n')
}
d<-do.call(rbind,ans);d$q_across_21_contrasts<-ave(d$pvalue,interaction(d$compartment,d$level),FUN=function(p)p.adjust(p,'BH'));write.csv(d,file.path(out,'taxonomy_pairwise_DESeq2.csv'),row.names=FALSE);write.csv(do.call(rbind,diag),file.path(out,'taxonomy_model_diagnostics_DESeq2.csv'),row.names=FALSE)
# Native R and PyDESeq2 need not produce identical dispersion estimates or P values.
pyfile<-file.path(root,'python_results','DA','taxonomy_pairwise_PyDESeq2.csv')
if(file.exists(pyfile)){p<-read.csv(pyfile,check.names=FALSE);m<-merge(p,d,by=c('compartment','level','popA','popB','feature'),suffixes=c('_Python','_R'));z<-lapply(split(m,interaction(m$compartment,m$level)),function(x){ok<-is.finite(x$log2FoldChange_Python)&is.finite(x$log2FoldChange_R);data.frame(compartment=x$compartment[1],level=x$level[1],n_coefficients=sum(ok),median_absolute_LFC_difference=median(abs(x$log2FoldChange_Python[ok]-x$log2FoldChange_R[ok])),n_significant_Python=sum(x$q_across_21_contrasts_Python<.05,na.rm=TRUE),n_significant_R=sum(x$q_across_21_contrasts_R<.05,na.rm=TRUE))});write.csv(do.call(rbind,z),file.path(out,'PyDESeq2_native_DESeq2_comparison.csv'),row.names=FALSE)}

 }), envir=new.env(parent=globalenv()))

 
 # Embedded module: protein_taxonomy.R
 eval(quote({
# Missing analysis: taxon-assigned protein NSAF, class/order/family.
suppressPackageStartupMessages(library(limma))
args<-commandArgs(trailingOnly=TRUE);root<-if(length(args))normalizePath(args[1]) else normalizePath('.')
out<-file.path(root,'R_results','DA');dir.create(out,recursive=TRUE,showWarnings=FALSE)
meta<-read.csv(file.path(root,'input','meta_data_TZ_AU.csv'),check.names=FALSE);prot<-read.csv(file.path(root,'input','metaproteome_TZ_AU.csv'),check.names=FALSE)
samples<-meta$id;pops<-unique(meta$group);g<-factor(meta$group,levels=pops);design<-model.matrix(~0+g);colnames(design)<-pops
prot$partner<-'Other';prot$partner[prot$superKingdom=='Bacteria' & prot$Phylum!='Cyanobacteria']<-'Noncyano_bacteria';prot$partner[prot$Phylum=='Ascomycota']<-'Mycobiont';prot$partner[prot$Phylum=='Chlorophyta']<-'Phycobiont';prot$partner[prot$superKingdom=='Bacteria' & prot$Phylum=='Cyanobacteria']<-'Cyanobiont'
tax<-read.csv(file.path(root,'input','amplicon_16S_tax_TZ_AU.csv'),row.names=1,check.names=FALSE);tax[is.na(tax)]<-'unclassified';asv<-as.matrix(read.csv(file.path(root,'input','amplicon_16S_asv_TZ_AU.csv'),row.names=1,check.names=FALSE)[,samples]);tax<-tax[rownames(asv),,drop=FALSE];org<-apply(tax,1,function(z)any(grepl('chloroplast|mitochond',z,ignore.case=TRUE)));core<-rowSums(asv>0)==21 & tax$Kingdom=='Bacteria' & !org;corefam<-setdiff(unique(tax$Family[core]),c('unclassified','various'))
masks<-list(Whole=rep(TRUE,nrow(prot)),All_bacteria=prot$superKingdom=='Bacteria',Noncyano_bacteria=prot$partner=='Noncyano_bacteria',Cyanobiont=prot$partner=='Cyanobiont',Phycobiont=prot$partner=='Phycobiont',Mycobiont=prot$partner=='Mycobiont',Shared_core_lineages=prot$partner=='Noncyano_bacteria' & prot$Family %in% corefam)
stopifnot(sum(masks$Shared_core_lineages)==259)
moderate<-function(fit){n<-nrow(fit$coefficients);tryCatch(eBayes(fit,trend=n>=30,robust=n>=20),error=function(e){warning(paste('Empirical-Bayes fallback to no trend/nonrobust:',conditionMessage(e)));eBayes(fit,trend=FALSE,robust=FALSE)})}
results<-list();overall<-list();ccres<-list();diag<-list();coverage<-list()
for(comp in names(masks))for(lev in c('Family','Order','Class')){
 v<-as.matrix(prot[masks[[comp]],samples]);storage.mode(v)<-'double';lab<-as.character(prot[[lev]][masks[[comp]]]);lab[is.na(lab)|lab=='']<-'unclassified'
 M<-rowsum(v,lab,reorder=FALSE);M<-sweep(M,2,colSums(v),'/')
 write.csv(data.frame(feature=rownames(M),M,check.names=FALSE),file.path(out,paste0('protein_taxonomy_matrix_',comp,'_',lev,'.csv')),row.names=FALSE)
 named<-!tolower(rownames(M)) %in% c('unclassified','various','unknown','na','nan')
 coverage[[length(coverage)+1]]<-data.frame(compartment=comp,level=lev,sample=samples,population=as.character(g),classified_fraction=colSums(M[named,,drop=FALSE]),n_named_taxa=sum(named),n_named_protein_groups=sum(!tolower(lab)%in%c('unclassified','various','unknown','na','nan')))
 M<-M[named,,drop=FALSE]
 if(!nrow(M)){diag[[length(diag)+1]]<-data.frame(compartment=comp,level=lev,n_named_features=0,n_base_features=0,n_global_eligible=0);next}
 Y<-log2(M);Y[!is.finite(Y)]<-NA;nd<-do.call(cbind,lapply(pops,function(z)rowSums(M[,g==z,drop=FALSE]>0)));colnames(nd)<-pops;rownames(nd)<-rownames(M)
 df<-rowSums(is.finite(Y))-rowSums(nd>0);basekeep<-rowSums(is.finite(Y))>=6 & df>=3
 if(!sum(basekeep)){diag[[length(diag)+1]]<-data.frame(compartment=comp,level=lev,n_named_features=nrow(M),n_base_features=0,n_global_eligible=0);next}
 fit<-lmFit(Y[basekeep,,drop=FALSE],design);idsbase<-rownames(fit$coefficients);glob<-apply(nd[idsbase,,drop=FALSE]>=2,1,all)
 if(sum(glob)){
  contrasts<-sapply(2:7,function(j){z<-rep(0,7);z[1]<- -1;z[j]<-1;z});rownames(contrasts)<-pops;fo<-moderate(contrasts.fit(fit[glob,],contrasts));tt<-topTable(fo,coef=1:6,number=Inf,sort.by='none')
  overall[[length(overall)+1]]<-data.frame(compartment=comp,level=lev,feature=rownames(tt),F=tt$F,pvalue=tt$P.Value,padj=tt$adj.P.Val)
 }
 for(ii in combn(7,2,simplify=FALSE)){
  a<-pops[ii[1]];b<-pops[ii[2]];na<-nd[,a];nb<-nd[,b];ids<-rownames(M)[basekeep & na>=2 & nb>=2]
  if(!length(ids))next
  ct<-matrix(0,7,1,dimnames=list(pops,paste0(b,'_minus_',a)));ct[a,1]<- -1;ct[b,1]<-1
  fc<-moderate(contrasts.fit(fit[match(ids,idsbase),],ct));tt<-topTable(fc,number=Inf,sort.by='none',confint=TRUE)
  results[[length(results)+1]]<-data.frame(compartment=comp,level=lev,popA=a,popB=b,feature=rownames(tt),log2FoldChange=tt$logFC,CI_L=tt$CI.L,CI_R=tt$CI.R,t=tt$t,pvalue=tt$P.Value,padj=tt$adj.P.Val,n_detected_A=na[rownames(tt)],n_detected_B=nb[rownames(tt)])
  idscc<-ids[na[ids]==3 & nb[ids]==3]
  if(length(idscc)){
   fc<-moderate(contrasts.fit(fit[match(idscc,idsbase),],ct));tt<-topTable(fc,number=Inf,sort.by='none')
   ccres[[length(ccres)+1]]<-data.frame(compartment=comp,level=lev,popA=a,popB=b,feature=rownames(tt),log2FoldChange=tt$logFC,pvalue=tt$P.Value,padj=tt$adj.P.Val)
  }
 }
 diag[[length(diag)+1]]<-data.frame(compartment=comp,level=lev,n_named_features=nrow(M),n_base_features=sum(basekeep),n_global_eligible=sum(glob));cat('Completed',comp,lev,'\n');flush.console()
}
d<-do.call(rbind,results);d$q_across_21_contrasts<-ave(d$pvalue,interaction(d$compartment,d$level),FUN=function(x)p.adjust(x,'BH'));write.csv(d,file.path(out,'protein_taxonomy_pairwise_limma.csv'),row.names=FALSE)
cc<-do.call(rbind,ccres);cc$q_across_21_contrasts<-ave(cc$pvalue,interaction(cc$compartment,cc$level),FUN=function(x)p.adjust(x,'BH'));write.csv(cc,file.path(out,'protein_taxonomy_pairwise_completecase.csv'),row.names=FALSE)
o<-do.call(rbind,overall);o$q_across_omnibus<-p.adjust(o$pvalue,'BH');write.csv(o,file.path(out,'protein_taxonomy_omnibus_F.csv'),row.names=FALSE)
write.csv(do.call(rbind,coverage),file.path(out,'protein_taxonomy_assignment_coverage.csv'),row.names=FALSE);write.csv(do.call(rbind,diag),file.path(out,'protein_taxonomy_model_diagnostics.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'));cat('Protein taxonomy finished\n')
writeLines(capture.output(warnings()),file.path(out,'protein_taxonomy_warnings.txt'))

 }), envir=new.env(parent=globalenv()))

 publication_figures(environment())
 pycoef<-file.path(root,'python_results','protein_OLS_validation_only.csv');if(file.exists(pycoef)){aa<-read.csv(pycoef,check.names=FALSE);bb<-read.csv(file.path(out,'DA','functional_pairwise_limma.csv'),check.names=FALSE);cc<-merge(aa,bb,by=c('compartment','level','popA','popB','feature'),suffixes=c('_Python','_R'));err<-max(abs(cc$log2FoldChange_Python-cc$log2FoldChange_R),na.rm=TRUE);write(data.frame(n_compared=nrow(cc),max_absolute_LFC_error=err,passed=err<1e-8),'Python_R_limma_coefficient_check')}
 cat('\nKEY RESULTS\n');print(centrows[centrows$method=='full_PCoA_Lingoes' & ((centrows$dataset=='16S' & centrows$level=='Family') | (centrows$dataset=='Function' & centrows$level=='og' & centrows$compartment %in% c('Whole','All_bacteria'))),]);print(ft);cat('\nFULL POPULATION EFFECTS\n');print(global);cat('\nMANTEL RESULTS\n');print(mm);cat('\n16S FAMILY AND OG CENTROID COMPARISONS\n');print(sites[sites$method=='full_PCoA_Lingoes' & ((sites$dataset=='16S' & sites$level=='Family') | (sites$dataset=='Function' & sites$level=='og' & sites$compartment %in% c('Whole','All_bacteria'))),]);for(co in c('Cyanobiont','Mycobiont','Noncyano_bacteria','Phycobiont'))for(le in c('subRole','og')){cat('\nTOP INFORMATIVE FUNCTIONS: ',co,' ',le,'\n',sep='');zz<-rankings[rankings$partner==co & rankings$level==le & rankings$informative_annotation,];print(head(zz[order(zz$rank_informative),c('feature','description','mean_fraction','n_thalli_detected','n_populations_detected_2of3')],10))};print(sessionInfo())
 writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'));saveRDS(list(meta=meta,profiles=profiles,centroids=cent,global=global,core=ca,family_tests=ft,rankings=rankings),file.path(out,'analysis_objects.rds'))
 cat('\nALL R ANALYSES COMPLETED SUCCESSFULLY\nFinished: ',format(Sys.time()),'\n',sep='');writeLines('ALL R ANALYSES COMPLETED SUCCESSFULLY',file.path(out,'COMPLETED.txt'))
}

render_existing <- function() {
 needed<-c('vegan','ggplot2','ggrepel','limma','statmod');missing<-needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)];if(length(missing))stop('Missing packages: ',paste(missing,collapse=', '),'. Install these CRAN packages first.')
args<-commandArgs(trailingOnly=TRUE);root<-if(length(args))normalizePath(args[1]) else normalizePath(getwd());out<-file.path(root,'R_results');path<-file.path(out,'analysis_objects.rds')
if(!file.exists(path))stop('analysis_objects.rds is missing. Complete scripts/run_R.R first. For an interrupted run, send the existing console_log.txt and ERROR.txt.')
obj<-readRDS(path);meta<-obj$meta;samples<-meta$id;pops<-c('Mweka1','Mweka2','Mweka3','Marangu','Maua','Umbwe','Styria');g<-factor(meta$group,levels=pops);profiles<-obj$profiles;cent<-obj$centroids;ca<-obj$core;ft<-obj$family_tests;rankings<-obj$rankings
if(!'informative_annotation'%in%names(rankings)){rankings$informative_annotation<-!rankings$description%in%c('unclassified','various','Function unknown','poorly characterized') & !grepl('unknown function|domain of unknown|DUF[0-9]',rankings$description,ignore.case=TRUE);rankings$rank_informative<-NA_real_;ii<-which(rankings$informative_annotation);rankings$rank_informative[ii]<-ave(-rankings$mean_fraction[ii],interaction(rankings$partner[ii],rankings$level[ii]),FUN=function(x)rank(x,ties.method='min'))}
cc<-read.csv(file.path(out,'shared_family_contributions.csv'),check.names=FALSE);coverage<-split(cc,cc$family)
C<-as.matrix(read.csv(file.path(root,'input','amplicon_16S_asv_TZ_AU.csv'),row.names=1,check.names=FALSE)[,samples]);tax<-read.csv(file.path(root,'input','amplicon_16S_tax_TZ_AU.csv'),row.names=1,check.names=FALSE);tax<-tax[rownames(C),,drop=FALSE];tax[is.na(tax)]<-'unclassified';org<-apply(tax,1,function(z)any(grepl('chloroplast|mitochond',z,ignore.case=TRUE)));keep<-tax$Kingdom=='Bacteria' & !org;C<-C[keep,,drop=FALSE];tax<-tax[keep,,drop=FALSE]
unknown<-c('unclassified','various','Function unknown','poorly characterized','unknown','NA','');named<-function(x)!x%in%unknown;relative<-function(x)sweep(x,1,rowSums(x),'/');write<-function(x,n)write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE)

prot<-read.csv(file.path(root,'input','metaproteome_TZ_AU.csv'),check.names=FALSE);prot[is.na(prot)]<-'unclassified';V<-as.matrix(prot[,samples]);storage.mode(V)<-'double'
prot$partner<-'Other';prot$partner[prot$superKingdom=='Bacteria' & prot$Phylum!='Cyanobacteria']<-'Noncyano_bacteria';prot$partner[prot$Phylum=='Ascomycota']<-'Mycobiont';prot$partner[prot$Phylum=='Chlorophyta']<-'Phycobiont';prot$partner[prot$superKingdom=='Bacteria' & prot$Phylum=='Cyanobacteria']<-'Cyanobiont'
mm<-read.csv(file.path(out,'Mantel.csv'),check.names=FALSE)
publication_figures(environment())
}
run_selected <- function() {
 args<-commandArgs(trailingOnly=TRUE);root<-if(length(args))normalizePath(args[1],mustWork=TRUE) else normalizePath(getwd(),mustWork=TRUE)
 if(!RUN_MODE %in% c('auto','full','figures'))stop('RUN_MODE must be auto, full, or figures.')
 required<-c('analysis_objects.rds','COMPLETED.txt','Mantel.csv','shared_family_contributions.csv','DA/taxonomy_pairwise_DESeq2.csv','DA/functional_pairwise_limma.csv')
 complete<-all(file.exists(file.path(root,'R_results',required)))
 if(RUN_MODE=='figures' && !complete)stop('Completed R analysis objects and native model tables are required for figures mode. Set RUN_MODE to full to regenerate.')
 if(RUN_MODE=='full' || !complete)run_full() else {message('Using completed native R results; exporting JPG figure suite.');render_existing()}
}
tryCatch(run_selected(),error=function(e){args<-commandArgs(trailingOnly=TRUE);base<-if(length(args))args[1] else getwd();dir.create(file.path(base,'R_results'),recursive=TRUE,showWarnings=FALSE);writeLines(c(conditionMessage(e),capture.output(sys.calls())),file.path(base,'R_results','FIGURE_WORKFLOW_ERROR.txt'));stop(e)})
