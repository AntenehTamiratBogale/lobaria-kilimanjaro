# Validation and limitations

The earlier native R workflow completed in the user's R environment. All 21 supplied Python-R comparison checks passed. A separate check matched 77,234 limma log-fold-change coefficients against independent detected-log means. Additional checks verified finite-value BH adjustments and protein detection eligibility. The comparison scope is recorded in the supplied reference CSVs.

These checks do not independently reproduce empirical-Bayes limma P values. Native DESeq2 and PyDESeq2 inference can differ; final reporting should consistently use native R model outputs. Partial-NA coefficient warnings arise from non-detection in some groups and require attention to contrast eligibility.

The later plotting revisions and added whole-metaproteome OG model excluding Cyanobacteria had not been executed here when this repository was prepared. Earlier grammar checks are not graphical or statistical execution checks. The included CI performs parsing and basic input validation, not the full permutation and model workflow.

Repository preparation preserves the scientific workflow and adds an option-based mode switch and entry-point scripts. It does not certify submission readiness. Review collection dates, laboratory settings, accession availability, author declarations, final graphics and manuscript values before tagging a publication release.
