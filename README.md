# TaxShare
This package implements a unique/shared taxonomic analysis to compare microbial communities between two groups of samples at a user-defined taxonomic level.

The analysis produces a summary table reporting the proportion of taxa that are unique to each group and those shared between groups, along with group-level abundance summaries. Statistical significance is assessed using a Wilcoxon rank-sum test, with optional false discovery rate (FDR) correction.

Results can be visualized as stacked bar plots illustrating the relative contribution of group-unique and shared taxa, facilitating intuitive comparison of community structure between groups.

<img width="1732" height="1272" alt="image" src="https://github.com/user-attachments/assets/ffb87013-dee5-4e3b-b621-3538de797092" />

https://pubmed.ncbi.nlm.nih.gov/39998294/

This R package performs a unique/shared taxonomic comparison between two groups of samples at a user-specified taxonomic level (e.g., phylum, family, genus).

OTUs (or ASVs) are aggregated to the selected taxonomic rank, and each taxon is classified as:

Unique to group A

Unique to group B

Shared between groups, with directionality based on relative abundance

The resulting table reports the proportion and abundance of taxa in each category. Group differences are evaluated using a Wilcoxon rank-sum test, with optional FDR correction for multiple testing.

The package includes visualization tools that generate stacked bar plots summarizing mean (or user-defined) population sizes across unique and shared taxa, enabling clear interpretation of taxonomic patterns between groups.

# Typical Workflow
1. Run "tax.shared.table" to make table
2. Run "tax.shared.plot" to plot the results

* multicomp is run in one step
1. Run "tax.shared.multicomp" to make the plot



