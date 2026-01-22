#' Row-wise Medians with Suppressed Warnings
#'
#' Computes the median of each row in a numeric matrix or data frame,
#' coercing values to numeric and suppressing warnings (e.g., from NAs or non-numeric values).
#'
#' @param m A matrix or data frame with numeric or coercible values.
#'
#' @return A numeric vector of row medians.
#' @examples
#' mat <- matrix(c(1, 2, NA, 4, 5, 6), nrow = 2)
#' rowMedians(mat)
#'
#' @export
rowMedians <- function(m) {
  apply(m, 1, function(x) suppressWarnings(median(as.numeric(x), na.rm = TRUE)))
}

#' Row-wise Means with Suppressed Warnings
#'
#' Computes the mean of each row in a numeric matrix or data frame,
#' coercing values to numeric and suppressing warnings.
#'
#' @param m A matrix or data frame with numeric or coercible values.
#'
#' @return A numeric vector of row means.
#' @examples
#' mat <- matrix(c(1, 2, NA, 4, 5, 6), nrow = 2)
#' rowMeans2(mat)
#'
#' @export
rowMeans2 <- function(m) {
  apply(m, 1, function(x) suppressWarnings(mean(as.numeric(x), na.rm = TRUE)))
}

#' Determines which taxa are Shared and Unique Taxa Between Groups
#'
#' Compares OTU (Operational Taxonomic Unit) tables between two groups to determine
#' shared and unique taxa at a specified taxonomic level, and computes summary statistics
#' including p-values from Wilcoxon tests. Supports presence/absence or abundance data,
#' log-transformation, and optional inclusion of specific taxa.
#'
#' @param otu.input A numeric matrix or data frame of OTU counts with samples as rows and taxa as columns.
#' @param var A factor or character vector indicating group membership for each sample (length equal to number of rows in `otu.input`).
#' @param group1 Name of the first group (as in `var`) to compare.
#' @param group2 Name of the second group (as in `var`) to compare.
#' @param taxa2include Optional character vector of taxa names to include in the final output (pattern matching allowed).
#' @param tax.level A vector mapping each column of `otu.input` to its taxonomic family (or desired taxonomic level).
#' @param PresAbs Logical; if `TRUE`, use presence/absence data (default `FALSE`).
#' @param log Logical; if `TRUE`, apply log10 transformation to output (default `TRUE`). Ignored with a warning if `PresAbs = TRUE`.
#' @param min.obs Integer; minimum number of observations to consider a taxon "present" (default = 1).
#' @param round.to Integer; number of decimal places to round p-values (default = 4).
#' @param fdr Logical; if `TRUE`, apply FDR correction to p-values (default `TRUE`).
#' @param use.median Logical; if `TRUE`, use medians for summary statistics, otherwise use means (default `TRUE`).
#' @param add.n Logical; if `TRUE`, annotate columns with the number of contributing OTUs (default `TRUE`).
#'
#' @return A matrix of summarized abundance values by group and taxon, annotated with statistical significance of shared/unique presence.
#'
#' @details
#' This function compares OTU profiles between two groups by decomposing taxa into:
#' - Unique to group 1
#' - Shared between groups (group 1 abundance)
#' - Shared between groups (group 2 abundance)
#' - Unique to group 2
#'
#' It applies Wilcoxon rank-sum tests to assess statistical significance of shared and unique taxa, and reports adjusted p-values if desired.
#'
#' If `taxa2include` is supplied, only those taxa (or matching patterns) are retained in the final output.
#'
#' @examples
#' \dontrun{
#'
#' tax.result <- tax.shared.table(otu_matrix, sample_metadata$Group,
#'                                "Healthy", "Diseased", taxa2include = c("Bacteroides", "Prevotella"),
#'                                tax.level = tax$Family)
#' }
#'
#' @export

tax.shared.table <- function (otu.input,
                              var,
                              group1,
                              group2,
                              taxa2include,
                              tax.level = tax$family,
                              PresAbs = F,
                              min.obs = 1,
                              round.to = 3,
                              fdr = T,
                              summary.fun = mean,
                              add.n = T) {

  var <- factor(var)
  group1.otus <- otu.input[var == group1, ]
  group2.otus <- otu.input[var == group2, ]

  if (PresAbs) {
    gg <- "number of discrete OTUs "
    group1.otus[group1.otus >= min.obs] <- 1
    group2.otus[group2.otus >= min.obs] <- 1
  }else{
    gg <- ""
  }
  group1.otus.shared <- group1.otus.unique <- group1.otus
  group2.otus.shared <- group2.otus.unique <- group2.otus

  group1.otus.unique[, colSums(group2.otus) >= min.obs] <- 0
  group1.otus.shared[, colSums(group1.otus.unique) > 0] <- 0
  group2.otus.unique[, colSums(group1.otus) >= min.obs] <- 0
  group2.otus.shared[, colSums(group2.otus.unique) > 0] <- 0

  group1.fam.unique <- aggregate(t(group1.otus.unique), list(tax.level), sum)
  group1.fam.shared <- aggregate(t(group1.otus.shared), list(tax.level), sum)
  group2.fam.shared <- aggregate(t(group2.otus.shared), list(tax.level), sum)
  group2.fam.unique <- aggregate(t(group2.otus.unique), list(tax.level), sum)

  #####
  rownames(group1.fam.unique) <- group1.fam.unique$Group.1
  group1.fam.unique$Group.1 <- NULL
  g1.u <- apply(group1.fam.unique, 1, summary.fun, na.rm = TRUE)

  rownames(group1.fam.shared) <- group1.fam.shared$Group.1
  group1.fam.shared$Group.1 <- NULL
  g1.s <- apply(group1.fam.shared, 1, summary.fun, na.rm = TRUE)

  rownames(group2.fam.shared) <- group2.fam.shared$Group.1
  group2.fam.shared$Group.1 <- NULL
  g2.s <- apply(group2.fam.shared, 1, summary.fun, na.rm = TRUE)

  rownames(group2.fam.unique) <- group2.fam.unique$Group.1
  group2.fam.unique$Group.1 <- NULL
  g2.u <- apply(group2.fam.unique, 1, summary.fun, na.rm = TRUE)

  ps.unique.l <- ps.unique.g <- ps.shared <- shared.ef <-  NULL

  for (i in 1:nrow(group1.fam.unique)) {

    ps.unique.l[i] <- wilcox.test(t(group2.fam.unique[i, -1]),
                                  rep(0, nrow(group1.otus)),
                                  exact = F)$p.value

    ps.shared[i] <- wilcox.test(t(group1.fam.shared[i, -1]),
                                t(group2.fam.shared[i, -1]),
                                exact = F)$p.value

    if(median(t(group1.fam.shared[i, -1])) > median(t(group2.fam.shared[i, -1]))){
      shared.ef[i] <- group1
    }else{
      shared.ef[i] <- group2
    }

    ps.unique.g[i] <- wilcox.test(t(group1.fam.unique[i, -1]),
                                  rep(0, nrow(group2.otus)),
                                  exact = F)$p.value
  }

  if(fdr){
    ps.unique.l <- p.adjust(ps.unique.l, method = "fdr", n = nrow(group1.fam.unique))
    ps.shared <- p.adjust(ps.shared, method = "fdr", n = nrow(group1.fam.unique))
    ps.unique.g <- p.adjust(ps.unique.g, method = "fdr",n = nrow(group1.fam.unique))
  }

  ps.unique.l <- round(ps.unique.l, round.to)
  ps.shared <- round(ps.shared, round.to)
  ps.unique.g <- round(ps.unique.g, round.to)

  hh <- which(ps.unique.g > 0.05 | is.na(ps.unique.g))
  ps.unique.g[hh] <- " "
  ps.unique.g[-hh] <- paste0("(unique to ", group1, " p=", ps.unique.g[-hh], ")")

  kk <- which(ps.shared > 0.05 | is.na(ps.shared))
  ps.shared[kk] <- " "
  ps.shared[-kk] <- paste0("(shared ", shared.ef[-kk], " favored p=", ps.shared[-kk],")")

  jj <- which(ps.unique.l > 0.05 | is.na(ps.unique.l))
  ps.unique.l[jj] <- " "
  ps.unique.l[-jj] <- paste0("(unique to ", group2, " p=", ps.unique.l[-jj], ")")


  tax.res <- rbind(g1.u, g1.s, g2.s, g2.u)

  if(add.n){
    #otu.pick <- otu.input[var %in% c(group1, group2), ]
    #ns <- table(tax.level[ colSums(otu.pick) > 0]) # number of taxa in both groups
    #colnames(tax.res)[match(names(ns), colnames(tax.res))] <- paste0(names(ns), "(nOTUs=", ns, ")")

    tax.level <- as.factor(tax.level)
    ns.1u <- table(tax.level[ colSums(group1.otus.unique) > 0])
    ns.1s <- table(tax.level[ colSums(group1.otus.shared) > 0])
    ns.2u <- table(tax.level[ colSums(group2.otus.unique) > 0])

    colnames(tax.res)[match(names(ns.1u), colnames(tax.res))] <- paste0(names(ns.1u), " & (nOTUs=", ns.1u, ", ", ns.1s,   ", ", ns.2u, ") && ")
  }

  colnames(tax.res) <- paste0(colnames(tax.res),  ps.unique.g,
                              ps.shared, ps.unique.l)

  rownames(tax.res) <- c(paste0("unique to ", group1),
                         paste0("shared - ", gg, "in ", group1),
                         paste0("shared - ", gg, "in ", group2),
                         paste0("unique to ", group2))

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

  res <- NULL
  res$table <- tax.res
  res$summary.fun <- deparse(substitute(summary.fun))
  res$fdr <- fdr

  return(res)
}


#' Visualize Shared and Unique Taxa Between Groups
#'
#' This function visualizes differences in taxonomic abundance between two groups
#' by creating a horizontal barplot of shared and unique taxa, based on the results from
#' `tax.shared.table()`. It supports log transformation and subsetting to specific taxa.
#'
#' @param tax.res result from tax.shared.table
#' @param log Logical; if TRUE, apply log10(x + 1) transformation to the matrix and plot axis. Default is TRUE.
#'
#' @details This function builds upon `tax.shared.table()` to generate a visual summary of taxa shared and unique to two groups. Bar lengths represent median abundance (or presence) for each taxon in the respective group.
#'
#' @importFrom graphics barplot legend axis segments text par
#' @export
tax.shared.plot <- function (tax.shared.table_result,
                             log = T,
                             min.obs = 1,
                             color1 = "orange" ,# colors of var1.groups group 1
                             color2 = "purple" ,# colors of var1.groups group 2
                             family.name.length = 2,
                             round.to =  3) {

  tax.res <- tax.shared.table_result$table

  max.bar.size <- max(tax.res)

  if (log) {
    tax.res <- log10(as.matrix(tax.res) + 1)
  }

  par(xpd = TRUE)
  tax.res <- tax.res[, order(colSums(tax.res), decreasing = F)]
  #here <- ceiling(max(colSums(tax.res))) -.5
  here <- max(colSums(tax.res))

  ##################
  shared.1 <- colorspace::mixcolor(0.3, colorspace::sRGB(t(col2rgb(color1))),
                                   colorspace::sRGB(t(col2rgb(color2))))


  shared.2 <- colorspace::mixcolor(0.7, colorspace::sRGB(t(col2rgb(color1))),
                                   colorspace::sRGB(t(col2rgb(color2))))

  shared.1.col <- rgb(shared.1@coords, maxColorValue = 255)
  shared.2.col <- rgb(shared.2@coords, maxColorValue = 255)

  a <- 1
  b <- 2

  while(all(a != b)){
    a <- colnames(tax.res) <- gsub(") ", ")", colnames(tax.res))
    b <- gsub(") ", ")", a)
  }

  labs <- gsub(")(", ") (", b, fixed = T)

  barplot(tax.res,
          cex.axis = 0.75,
          cex.names = 0.75,
          font.lab = 4,
          xaxt = "n",
          yaxt = "n",
          horiz = T,
          main = "",
          xlab = paste(tax.shared.table_result$summary.fun ,"population size",
                       "\nFDR=",tax.shared.table_result$fdr),
          las = 1,
          col = c(color1, shared.1.col,
                  shared.2.col, color2))



  if(family.name.length == 1){
    labs <- sapply(strsplit(labs, " & "), "[[", 1)
  }

  if(family.name.length == 2){
    labs <- sapply(strsplit(labs, "&& "), "[[", 1)
    labs <- gsub("&  ", "", labs)
  }
  if(family.name.length == 3){

    labs <- gsub("&", "", labs)
    labs <- gsub("&&", "", labs)
    labs <- gsub("  ", " ", labs)
    labs <- gsub("  ", " ", labs)
  }

  axis(2,
       at= (1:length(labs) * 1.2) - .5,
       labels=labs,
       las = 2,
       cex.axis=.8,
       font = 4)

  usr <- par("usr")

  ytop <- mean(c(usr[1], usr[4]))
  y.adj <- diff(c(usr[1], ytop)) / 15

  counts <- c(1, 5, 10, 50, 100, 1000)
  labels <- paste0(counts, "\n")
  labels[1] <- paste0("Counts\n", labels[1], "\n")

  labels <- labels[counts < max.bar.size]
  counts <- counts[counts < max.bar.size]

  for (i in seq_along(counts)) {

    y <- ytop - (i - 1) * y.adj

    if (log ) {

      segments( x0 = here - log10(counts[i] + 1),
                x1 = here,
                y0 = y,
                y1 = y,
                lwd = 2 )
    }else{
      segments( x0 = here - counts[i],
                x1 = here,
                y0 = y,
                y1 = y,
                lwd = 2 )
    }

    text(here,
         y,
         labels[i],
         cex = 0.5,
         font = 4)
  }


  legend("bottomright",
         fill = c(color1, shared.1.col,
                  shared.2.col, color2),
         legend = row.names(tax.res),
         bty = "n",
         text.font = 4,
         cex = 0.75)
}

tax.shared.multicomp <- function (otu.input ,
                                  var1 ,
                                  var2 ,
                                  var1.groups ,
                                  var2.groups ,
                                  tax.level ,
                                  taxa2include,
                                  family.name.length = 3,
                                  var1.color1 = "orange" ,# colors of var1.groups group 1
                                  var1.color2 = "purple" ,# colors of var1.groups group 2
                                  PresAbs = F ,
                                  log = T ,
                                  min.obs = 1,
                                  fdr = T,
                                  summary.fun = mean){

  if (length(var1.color1) == 1) {
    var1.color1 <- rep(var1.color1, nrow(otu.input))
  }
  if (length(var1.color2) == 1) {
    var1.color2 <- rep(var1.color2, nrow(otu.input))
  }

  var1.levels <- factor(var1.groups)
  if(length(var1.levels) != 2){
    stop("must compare EXACTLY 2 groups in var 1")
  }

  n <- length(var2.groups)

  if (missing(taxa2include)) {
    n.taxa <- length(unique(tax.level))
    taxa2include <- unique(tax.level)
  }else{
    n.taxa <- length(taxa2include)
  }

  add <- vector(length = n)
  add[-1] <- T
  res.all <- list()
  n <- length(var2.groups)

  for (i in 1:n) {
    print(i)

    var1.p <- var1[var2 == var2.groups[i]]
    otu.p <- otu.input[var2 == var2.groups[i], ]

    var1.p.groups <- unique(var1.p)

    if(length(unique(var1.p)) != 2){
      warning(var2.groups[i], " Does Not Contain 2 levels")
      res.all[[i]] <- NA

    }else{
      res.i <- tax.shared.table(otu.input = otu.p ,
                                var =  var1.p ,
                                group1 = var1.levels[1] ,
                                group2 = var1.levels[2] ,
                                tax.level = tax.level,
                                taxa2include = taxa2include,
                                PresAbs = PresAbs,
                                round.to = 3,
                                add.n = T,
                                min.obs = min.obs,
                                fdr = fdr,
                                summary.fun = summary.fun)$table

      if (log) {
        res.i <- log10(as.matrix(res.i) + 1)
      }

      labs <- colnames(res.i)

      if(family.name.length == 1){
        labs <- sapply(strsplit(labs, " & "), "[[", 1)
      }

      if(family.name.length == 2){
        labs <- sapply(strsplit(labs, "&& "), "[[", 1)
        labs <- gsub("&  ", "", labs)
      }
      if(family.name.length == 3){

        labs <- gsub("&", "", labs)
        labs <- gsub("  ", " ", labs)
        labs <- gsub("  ", " ", labs)
      }

      colnames(res.i) <- labs
      ########
      res.all[[i]] <- res.i
    }
  }

  if (log) {
    max.bar.size <- 10 ^ max(unlist(res.all))
  }else{
    max.bar.size <- max(unlist(res.all))
  }

  names(res.all) <- var2.groups

  complete.comps <- !is.na(res.all)
  #  res.all <- res.all[complete.comps]

  here = max(sapply(res.all[complete.comps], function(x) colSums(x, na.rm = T)))
  here <- here * 1.5
  unique.1.col <- NULL
  shared.1.col <- NULL
  shared.2.col <- NULL
  unique.2.col <- NULL

  par(xpd = TRUE)

  for (i in (1 : n)[complete.comps]) {
    p <- var2 == var2.groups[i]
    unique.1.col[i] <- var1.color1[p][1]
    unique.2.col[i] <- var1.color2[p][1]
    if (i != 1) {
      new.colnames <- gsub(paste0(taxa2include, collapse = "|"),
                           "", colnames(res.all[[i]]))
      new.colnames <- gsub("  ", "", new.colnames)
    }
    else {
      new.colnames <- colnames(res.all[[i]])
    }
    shared.1 <- colorspace::mixcolor(0.3, colorspace::sRGB(t(col2rgb(var1.color1[p][1]))),
                                     colorspace::sRGB(t(col2rgb(var1.color2[p][1]))))
    shared.2 <- colorspace::mixcolor(0.7, colorspace::sRGB(t(col2rgb(var1.color1[p][1]))),
                                     colorspace::sRGB(t(col2rgb(var1.color2[p][1]))))
    shared.1.col[i] <- rgb(shared.1@coords, maxColorValue = 255)
    shared.2.col[i] <- rgb(shared.2@coords, maxColorValue = 255)
    colnames(res.all[[i]]) <- new.colnames

    x.lab <- paste(deparse(substitute(summary.fun)), "population size\nFDR =", fdr)

    barplot(res.all[[i]], ylim = c(0, ncol(res.all[[i]]) *
                                     (n + 1)), xlim = c(0, here),
            cex.axis = 0.5, cex.names = 0.5,
            xaxt = "n", horiz = T, main = "",
            font.lab = 4,
            xlab = x.lab,
            las = 1, add = add[i], space = c(n - i, rep(n, ncol(res.all[[i]]))[-1]),
            col = c(unique.1.col[i],
                    shared.1.col[i],
                    shared.2.col[i],
                    unique.2.col[i]))
  }

  usr <- par("usr")

  ytop <- mean(c(usr[1], usr[4]))
  y.adj <- diff(c(usr[1], ytop)) / 15

  counts <- c(1, 5, 10, 50, 100, 1000)
  labels <- paste0(counts, "\n")
  labels[1] <- paste0("Counts\n", labels[1], "\n")

  labels <- labels[counts < max.bar.size]
  counts <- counts[counts < max.bar.size]

  for (i in seq_along(counts)) {

    y <- ytop - (i - 1) * y.adj

    if (log ) {

      segments( x0 = here - log10(counts[i] + 1),
                x1 = here,
                y0 = y,
                y1 = y,
                lwd = 2 )
    }else{
      segments( x0 = here - counts[i],
                x1 = here,
                y0 = y,
                y1 = y,
                lwd = 2 )
    }

    text(here,
         y,
         labels[i],
         cex = 0.5,
         font = 4)
  }

  add.legend(here, n.taxa,
             unique.1.col, shared.1.col, shared.2.col, unique.2.col,
             var1.groups, var2.groups)

  return(res.all)
}

add.legend <- function(here, n.taxa,
                       unique.1.col, shared.1.col, shared.2.col, unique.2.col,
                       var1.groups, var2.groups){

  all_colors <- c(unique.1.col, shared.1.col, shared.2.col, unique.2.col)

  # Position parameters
  x0 <- here + here/5  # left edge
  y0 <- n.taxa # top edge
  dx <- 2.5   # box width
  dy <- .5   # box height
  ncol <- 4   # 4 columns
  nrow <- length(unique.1.col) # number of groups in 2nd category

  # Draw boxes
  for (i in 1 : length(all_colors)) {
    row <- (i - 1) %% nrow
    col <- (i - 1) %/% nrow

    xleft <- x0 + col * dx
    ytop  <- y0 - row * dy

    rect(xleft, ytop - dy, xleft + dx, ytop,
         col = all_colors[i], border = "black")

  }

  x.leg.lables <- c("U1", "S1" , "S2", "U2")

  for (i in 1 : 4) {

    col <- (0:3)[i]

    xleft <- x0 + col * dx

    text(xleft + dx * 0.5, y0 + .25,
         labels = x.leg.lables[i], cex = 0.5)
  }

  for (i in 1 : length(var2.groups)) {

    row <- (0:length(var2.groups))[i]

    xleft <- x0 + 4 * dx
    ytop  <- y0 - row * dy
    text(xleft, ytop - dy * 0.5,
         labels = var2.groups[i], cex = 0.5, pos = 4)
  }

  x.leg.lables2 <- c(var1.groups[1],
                     "---------",
                     var1.groups[2])
  for (i in 1:3){

    col <- c(0,1.5,3)[i]

    xleft <- x0 + col * dx

    print(xleft + dx * 0.5)

    text(xleft + dx * 0.5, y0 - (length(var2.groups) * dy) - .25,
         labels = x.leg.lables2[i], cex = .5)
  }
}



