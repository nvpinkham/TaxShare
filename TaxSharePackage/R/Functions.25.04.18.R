#' Compare Shared and Unique Taxa Between Groups
#'
#' This function compares OTU abundances between two groups
#' and summarizes shared and unique taxa at a given taxonomic level. It can apply presence-absence
#' filtering and log transformation, and optionally include p-values for group differences.
#'
#' @param otu A matrix or data frame of OTU abundances (samples in rows, OTUs in columns).
#' @param var A vector indicating group membership for each sample (must match rows of `otu`).
#' @param group1 A value from `var` representing the first group.
#' @param group2 A value from `var` representing the second group.
#' @param tax.level A vector assigning taxonomic classification (e.g., family) to each OTU column.
#' @param PresAbs Logical; if TRUE, convert abundances to presence/absence (binary) data. Default is FALSE.
#' @param log Logical; if TRUE, apply log10(x + 1) transformation to final output. Default is TRUE.
#' @param taxa2include Optional vector of taxa names to include in the final output. If provided, will subset and order columns.
#' @param round.to Number of decimal places for rounding adjusted p-values. Default is 3.
#' @param add.n Logical; if TRUE, adds the number of OTUs per taxon to column names. Default is TRUE.
#'
#' @return A matrix with 4 rows (unique and shared taxa for each group) and columns representing taxa,
#' possibly with attached p-values and OTU counts. The matrix is log-transformed if `log = TRUE`.
#'
#' @details The function aggregates OTU counts to the specified taxonomic level, identifies taxa unique to each group,
#' and computes Wilcoxon tests for statistical differences in abundance. Shared and unique taxa are tested separately,
#' and FDR-adjusted p-values are included in column labels.
#'
#' @importFrom stats wilcox.test p.adjust
#' @importFrom matrixStats rowMedians
#' @export
tax.shared.table <- function (otu, var, group1, group2, tax.level, PresAbs = F,
                              log = T, taxa2include, round.to = 3, add.n = T)
{
  if (PresAbs & log) {
    warning("Log transformation not informative when using PresAbs")
  }
  var <- as.factor(var)
  group1.otus <- otu[var == group1, ]
  group2.otus <- otu[var == group2, ]
  if (PresAbs) {
    gg <- "number of discrete OTUs"
    group1.otus[group1.otus > 0] <- 1
    group2.otus[group2.otus > 0] <- 1
  }else {
    gg <- ""
  }
  group1.otus.shared <- group1.otus.unique <- group1.otus
  group2.otus.shared <- group2.otus.unique <- group2.otus
  group1.otus.unique[, colSums(group2.otus) > 0] <- 0
  group1.otus.shared[, colSums(group1.otus.unique) > 0] <- 0
  group2.otus.unique[, colSums(group1.otus) > 0] <- 0
  group2.otus.shared[, colSums(group2.otus.unique) > 0] <- 0
  group1.fam.unique <- aggregate(t(group1.otus.unique), list(tax.level),
                                 sum)
  group1.fam.shared <- aggregate(t(group1.otus.shared), list(tax.level),
                                 sum)
  group2.fam.unique <- aggregate(t(group2.otus.unique), list(tax.level),
                                 sum)
  group2.fam.shared <- aggregate(t(group2.otus.shared), list(tax.level),
                                 sum)
  g1.u <- rowMedians(group1.fam.unique[, -1])
  names(g1.u) <- group1.fam.unique$Group.1
  g1.s <- rowMedians(group1.fam.shared[, -1])
  names(g1.s) <- group1.fam.shared$Group.1
  g2.u <- rowMedians(group2.fam.unique[, -1])
  names(g2.u) <- group2.fam.unique$Group.1
  g2.s <- rowMedians(group2.fam.shared[, -1])
  names(g2.s) <- group2.fam.shared$Group.1
  ps.unique.l <- ps.unique.g <- ps.shared <- NULL
  for (i in 1:nrow(group1.fam.unique)) {
    ps.unique.l[i] <- wilcox.test(t(group1.fam.unique[i,
                                                      -1]), t(group2.fam.unique[i, -1]), alternative = "less",
                                  exact = F)$p.value
    ps.shared[i] <- wilcox.test(t(group1.fam.shared[i, -1]),
                                t(group2.fam.shared[i, -1]), exact = F)$p.value
    ps.unique.g[i] <- wilcox.test(t(group1.fam.unique[i,
                                                      -1]), t(group2.fam.unique[i, -1]), alternative = "greater",
                                  exact = F)$p.value
  }
  ps.unique.l <- p.adjust(ps.unique.l, method = "fdr")
  ps.shared <- p.adjust(ps.shared, method = "fdr")
  ps.unique.g <- p.adjust(ps.unique.g, method = "fdr")
  ps.unique.l <- round(ps.unique.l, round.to)
  ps.shared <- round(ps.shared, round.to)
  ps.unique.g <- round(ps.unique.g, round.to)
  hh <- which(ps.unique.g > 0.05 | is.na(ps.unique.g))
  ps.unique.g[hh] <- " "
  ps.unique.g[-hh] <- paste0("(unique greater p=", ps.unique.g[-hh],
                             ")")
  kk <- which(ps.shared > 0.05 | is.na(ps.shared))
  ps.shared[kk] <- " "
  ps.shared[-kk] <- paste0("(shared p=", ps.shared[-kk], ")")
  jj <- which(ps.unique.l > 0.05 | is.na(ps.unique.l))
  ps.unique.l[jj] <- " "
  ps.unique.l[-jj] <- paste0("(unique less p=", ps.unique.l[-jj],
                             ")")
  tax.res <- rbind(g1.u, g1.s, g2.s, g2.u)
  if (add.n) {
    otu.pick <- otu[var %in% c(group1, group2), ]
    ns <- table(tax.level[colSums(otu.pick) > 0])
    colnames(tax.res)[match(names(ns), colnames(tax.res))] <- paste0(names(ns),
                                                                     "(nOTUs=", ns, ")")
  }
  colnames(tax.res) <- paste0(colnames(tax.res), ps.unique.l,
                              ps.shared, ps.unique.g)
  rownames(tax.res) <- c(paste0("unique to ", group1), paste0("shared - ",
                                                             gg, "in", group1), paste0("shared - ", gg, "in", group2),
                         paste0("unique to ", group2))
  tax.res <- tax.res[, colSums(tax.res) != 0]
  if (log) {
    tax.res <- log10(as.matrix(tax.res) + 1)
  }
  tax.res <- tax.res[, order(colSums(tax.res), decreasing = F)]
  if (!missing(taxa2include)) {
    print(taxa2include)
    tax.res2 <- as.data.frame(matrix(nrow = 4, ncol = length(taxa2include)))
    row.names(tax.res2) <- row.names(tax.res)
    for (i in 1:length(taxa2include)) {
      aa <- grep(taxa2include[i], colnames(tax.res))
      if (length(aa) > 0) {
        tax.res2[, i] <- tax.res[, aa]
        colnames(tax.res2)[i] <- colnames(tax.res)[aa]
      }
      else {
        colnames(tax.res2)[i] <- taxa2include[i]
        tax.res2[, i] <- 0
      }
    }
    tax.res <- as.matrix(tax.res2)
  }
  colnames(tax.res) <- gsub("(", " (", colnames(tax.res), fixed = T)
  return(tax.res)
}

#' Visualize Shared and Unique Taxa Between Groups
#'
#' This function visualizes differences in taxonomic abundance between two groups
#' by creating a horizontal barplot of shared and unique taxa, based on the results from
#' `tax.shared.table()`. It supports log transformation and subsetting to specific taxa.
#'
#' @param otu A matrix or data frame of OTU abundances (samples in rows, OTUs in columns).
#' @param var A vector indicating group membership for each sample (must match rows of `otu`).
#' @param group1 A value from `var` representing the first group.
#' @param group2 A value from `var` representing the second group.
#' @param tax.level A vector assigning taxonomic classification (e.g., family) to each OTU column. Default is `tax$family`.
#' @param PresAbs Logical; if TRUE, convert OTU counts to presence/absence (binary). Default is FALSE.
#' @param log Logical; if TRUE, apply log10(x + 1) transformation to the matrix and plot axis. Default is TRUE.
#' @param taxa2include Optional vector of taxa names to display in the plot. If supplied, only these taxa will be shown.
#'
#' @return A matrix of taxonomic group medians (possibly log-transformed), invisibly returned. Also generates a barplot for visual inspection.
#'
#' @details This function builds upon `tax.shared.table()` to generate a visual summary of taxa shared and unique to two groups. Bar lengths represent median abundance (or presence) for each taxon in the respective group.
#'
#' @importFrom graphics barplot legend axis segments text par
#' @export
tax.shared <- function (otu, var, group1, group2, tax.level = tax$family, PresAbs = F,
                        log = T, taxa2include){
  tax.res <- tax.shared.table(otu, var, group1, group2, tax.level,
                              PresAbs, log, taxa2include)
  par(mar = c(3, 20, 4, 4))
  par(xpd = TRUE)
  tax.res <- tax.res[, order(colSums(tax.res), decreasing = F)]
  here <- ceiling(max(colSums(tax.res))) - 0.5
  if (!missing(taxa2include)) {
    print(taxa2include)
    tax.res2 <- as.data.frame(matrix(nrow = 4, ncol = length(taxa2include)))
    row.names(tax.res2) <- row.names(tax.res)
    for (i in 1:length(taxa2include)) {
      aa <- grep(taxa2include[i], colnames(tax.res))
      if (length(aa) > 0) {
        tax.res2[, i] <- tax.res[, aa]
        colnames(tax.res2)[i] <- colnames(tax.res)[aa]
      }
      else {
        colnames(tax.res2)[i] <- taxa2include[i]
        tax.res2[, i] <- 0
      }
    }
    tax.res <- as.matrix(tax.res2)
  }
  a <- colnames(tax.res)
  b <- paste(colnames(tax.res), "dummy")
  while (all(a != b)) {
    a <- colnames(tax.res) <- gsub(") ", ")", colnames(tax.res))
    b <- gsub(") ", ")", a)
  }
  labs <- gsub(")(", ") (", b, fixed = T)
  barplot(tax.res, cex.axis = 0.75, cex.names = 0.75, font.lab = 4,
          xaxt = "n", yaxt = "n", horiz = T, main = "", xlab = "median population size (log 10)",
          las = 1, col = c("orange", "#b38a61", "#887382", "purple"))
  axis(2, at = (1:length(labs) * 1.2) - 0.5, labels = labs,
       las = 2, cex.axis = 0.5, font = 4)
  if (log) {
    segments(x0 = here - log10(10), x1 = here, y0 = 14, y1 = 14,
             lwd = 2)
    text(here - 0.5, 14, "Counts\n10\n\n", cex = 0.5, font = 4)
    segments(x0 = here - log10(100), x1 = here, y0 = 13,
             y1 = 13, lwd = 2)
    text(here - 0.5, 13, "100\n", cex = 0.5, font = 4)
    segments(x0 = here - 3, x1 = here, y0 = 12, y1 = 12,
             lwd = 2)
    text(here - 0.5, 12, "1000\n", cex = 0.5, font = 4)
  }
  legend("bottomright", fill = c("orange", "#b38a61", "#887382",
                                 "purple"), legend = row.names(tax.res), bty = "n", text.font = 4,
         cex = 0.75)
  return(tax.res)
}
