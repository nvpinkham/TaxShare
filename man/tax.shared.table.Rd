% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.25.06.17.R
\name{tax.shared.table}
\alias{tax.shared.table}
\title{Determines which taxa are Shared and Unique Taxa Between Groups}
\usage{
tax.shared.table(
  otu.input,
  var,
  group1,
  group2,
  taxa2include,
  tax.level = tax$family,
  PresAbs = F,
  log = T,
  min.obs = 1,
  round.to = 4,
  fdr = T,
  use.median = T,
  add.n = T
)
}
\arguments{
\item{otu.input}{A numeric matrix or data frame of OTU counts with samples as rows and taxa as columns.}

\item{var}{A factor or character vector indicating group membership for each sample (length equal to number of rows in \code{otu.input}).}

\item{group1}{Name of the first group (as in \code{var}) to compare.}

\item{group2}{Name of the second group (as in \code{var}) to compare.}

\item{taxa2include}{Optional character vector of taxa names to include in the final output (pattern matching allowed).}

\item{tax.level}{A vector mapping each column of \code{otu.input} to its taxonomic family (or desired taxonomic level).}

\item{PresAbs}{Logical; if \code{TRUE}, use presence/absence data (default \code{FALSE}).}

\item{log}{Logical; if \code{TRUE}, apply log10 transformation to output (default \code{TRUE}). Ignored with a warning if \code{PresAbs = TRUE}.}

\item{min.obs}{Integer; minimum number of observations to consider a taxon "present" (default = 1).}

\item{round.to}{Integer; number of decimal places to round p-values (default = 4).}

\item{fdr}{Logical; if \code{TRUE}, apply FDR correction to p-values (default \code{TRUE}).}

\item{use.median}{Logical; if \code{TRUE}, use medians for summary statistics, otherwise use means (default \code{TRUE}).}

\item{add.n}{Logical; if \code{TRUE}, annotate columns with the number of contributing OTUs (default \code{TRUE}).}
}
\value{
A matrix of summarized abundance values by group and taxon, annotated with statistical significance of shared/unique presence.
}
\description{
Compares OTU (Operational Taxonomic Unit) tables between two groups to determine
shared and unique taxa at a specified taxonomic level, and computes summary statistics
including p-values from Wilcoxon tests. Supports presence/absence or abundance data,
log-transformation, and optional inclusion of specific taxa.
}
\details{
This function compares OTU profiles between two groups by decomposing taxa into:
\itemize{
\item Unique to group 1
\item Shared between groups (group 1 abundance)
\item Shared between groups (group 2 abundance)
\item Unique to group 2
}

It applies Wilcoxon rank-sum tests to assess statistical significance of shared and unique taxa, and reports adjusted p-values if desired.

If \code{taxa2include} is supplied, only those taxa (or matching patterns) are retained in the final output.
}
\examples{
\dontrun{

tax.result <- tax.shared.table(otu_matrix, sample_metadata$Group,
                               "Healthy", "Diseased", taxa2include = c("Bacteroides", "Prevotella"),
                               tax.level = tax$Family)
}

}
