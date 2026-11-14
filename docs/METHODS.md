# Analysis definitions

## Core membership

Clean bacterial ASVs exclude organellar assignments. Primary core membership requires at least one read in every one of 18 Kilimanjaro thalli, independently of Styria. Recovery is then evaluated in Styria. Counts of 118, 92 and 109 refer respectively to the primary core, recovery in all three Styria thalli, and recovery in at least one Styria thallus in the previously completed analysis. Figures emphasize the first two sets. Threshold and rarefaction sensitivities remain available.

## Population structure

The scripts retain each of the seven populations separately. Bray-Curtis dissimilarities drive population tests and ordinations. Full analytical centroid distances use all positive PCoA dimensions after Lingoes correction; two-dimensional NMDS is visualization only. The seven possible exact reference choices imply a minimum reference-relabelling P of 1/7 = 0.143.

## Native feature models

16S taxonomic contrasts use raw aggregated counts with native DESeq2. Protein functional models use detected log2 relative NSAF in native limma, without imputation. Eligible protein contrasts require detection in at least two thalli in each compared population and sufficient information for model fitting. BH correction retains eligible features and all 21 population pairs within each compartment and resolution. The added whole-metaproteome OG exclusion sensitivity applies the same model rules after removing cyanobacterial protein groups.

## Additional analyses

Mantel tests use matrix permutations and explicitly defined scopes. Within-family functional analyses require adequate positive protein profiles; nonsignificance is not equivalence. Family-level alpha diversity aggregates named families, reports excluded unresolved bins, and uses one biological-thallus value after repeated common-depth rarefaction. Repeated draws do not increase the biological sample size.

## Graphics

The current workflow exports JPG only, with matching numbered captions. Taxonomy uses a shared categorical palette; NMDS uses filled circles without connecting lines; volcano displays use three panels and data-selected population contrasts. See `docs/FIGURES.md`. Rendered figures require final visual review.
