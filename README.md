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

Multliple comparisons is run in one step:
1. Run "tax.shared.multicomp" to make the plot

<div class="code-box">
  <div class="code-title">R Example</div>
  <pre><code>



map <- NULL
map$groups <- c("A", "A", "A", "A", "A", "A", "A", "A",
                "B", "B", "B", "B", "B",
                "C", "C", "C", "C")

map$season <- c("Spring", "Spring", "Spring","Spring", "Fall", "Fall", "Fall", "Fall",
                "Spring", "Spring", "Fall", "Fall", "Fall",
                "Spring", "Spring", "Fall", "Fall")

site <- as.data.frame(matrix(nrow = 17, ncol = 5))
colnames(site) <- c("Mule", "White_tail", "Black_Bear", "Beaver",
                    "Golden_Eagle")

rownames(site) <- paste0(rownames(site), ".", map$groups)

site$Mule <- c(25, 25, 14, 14, 0,0, 0, 0,
               0, 0, 0, 0, 0,
               21, 0, 11, 15)
#site$White_tail <- c(5, 5, 5, 5, 5, 5, 5, 5)
site$White_tail <- c(0, 0, 0, 0, 5, 5, 5, 4,
                     0, 0, 10, 5, 5,
                     0, 0, 4, 1)

site$Black_Bear <- c(1, 2, 1, 0, 1, 0,1, 0,
                     2, 1, 0, 0, 0,
                     1, 2, 0, 0)

site$Beaver <- c(0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 1, 1, 1,
                 0, 0, 0, 0)

#site$Golden_Eagle <- c(0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0)
#site$Golden_Eagle <- c(0, 1, 0, 0, 2, 1, 0, 0,  2, 3, 0, 0)
site$Golden_Eagle <- c(0, 0, 0, 0, 0, 0, 0, 0,
                       1, 1, 2, 0, 0,
                       0, 0, 1, 0)

#site$Golden_Eagle <- 1

tax <- NULL
tax$family <- c("Deer", "Deer", "Bear", "Beaver",
                "Eagle")
tax$class <- c("Mammalia", "Mammalia", "Mammalia", "Mammalia",
               "Aves")
tax$species <- colnames(site)

taxa2include <- unique(tax$family )

site[map$groups %in% c("A", "B" ), ]


tax.table <- tax.shared.table(otu.input = site,
                              var = map$groups,
                              group1 = "A",
                              group2 = "B",
                              min.obs = 1,
                              tax.level = tax$class,
                              PresAbs = F,
                              fdr = T,
                              round.to = 3,
                              summary.fun = mean)

par(mar = c(5, 10, 4, 4))
res <- tax.shared.plot(tax.shared.table_result = tax.table,
                       color1 = "salmon1",
                       color2 = "lightgreen",
                       log = F,
                       round.to = 3)

tax.table2 <- tax.shared.table(otu.input = site,
                               var = map$groups,
                               group1 = "A",
                               group2 = "B",
                               min.obs = 1,
                               tax.level = tax$family,
                               PresAbs = F,
                               fdr = F,
                               round.to = 3)

res2 <- tax.shared.plot(tax.shared.table_result = tax.table2,
                        color1 = "salmon1",
                        color2 = "lightgreen",
                        log = F,
                        round.to = 3)
    
  </code></pre>
</div>


---

---


