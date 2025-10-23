
library(TaxShare)

map <- NULL
map$groups <- c("A", "A", "A", "A", "A",
                "B", "B", "B", "B", "B")

site <- as.data.frame(matrix(nrow = 10, ncol = 4))
colnames(site) <- c("Mule", "White_tail", "Black_Bear", "Golden_Eagle")

site$Mule <- c(5, 5, 4, 5, 3,
               0, 0, 0, 0, 0)

site$White_tail <- c(3, 5, 6, 5, 1,
                     5, 0, 1, 4, 0)

site$Black_Bear <- c(1, 0, 0, 1, 0,
                     0, 2, 1, 0, 1)

site$Golden_Eagle <- c(0, 0, 0, 0, 0,
                       0, 1, 1, 2, 1)

tax <- NULL
tax$family <- c("Deer", "Deer", "Bear", "Eagle")
tax$species <- colnames(site)

tax.shared.plot(otu.input = site,
                var = map$groups,
                group1 = "A",
                group2 = "B",
                min.obs = 1,
                tax.level = tax$family,
                fdr = F,
                round.to = 3,
                log = F)


map <- NULL
map$groups <- c("A", "A", "A", "A", "A", "A", "A", "A",
                "B", "B", "B", "B", "B",
                "C", "C", "C", "C")

map$season <- c("Spring", "Spring", "Spring","Spring", "Fall", "Fall", "Fall", "Fall",
                "Spring", "Spring", "Fall", "Fall", "Fall",
                "Spring", "Spring", "Fall", "Fall")

site <- as.data.frame(matrix(nrow = 17, ncol = 4))
colnames(site) <- c("Mule", "White_tail", "Black_Bear", "Golden_Eagle")

site$Mule <- c(5, 5, 4, 4, 0,0, 0, 0,
               0, 0, 0, 0, 0,
               2, 0, 1, 5)

site$White_tail <- c(0, 0, 0, 0, 5, 5, 5, 4,
                     0, 0, 10, 5, 5,
                     0, 0, 4, 1)

site$Black_Bear <- c(1, 2, 1, 0, 1, 0,1, 0,
                     2, 1, 0, 0, 0,
                     1, 2, 0, 0)

site$Golden_Eagle <- c(0, 0, 0, 0, 0, 0, 0, 0,
                       1, 1, 2, 0, 0,
                       0, 0, 1, 0)

tax <- NULL
tax$family <- c("Deer", "Deer", "Bear", "Eagle")
tax$species <- colnames(site)

res <- tax.shared.plot(otu.input = site,
                       var = map$season,
                       group1 = "Fall",
                       group2 = "Spring",
                       min.obs = 1,
                       tax.level = tax$family,
                       PresAbs = F,
                       fdr = F,
                       round.to = 3)


p <- map$groups == "A"
res <- tax.shared.plot(otu.input = site[p, ],
                       var = map$season[p],
                       group1 = "Spring",
                       group2 = "Fall",
                       tax.level = tax$family,
                       PresAbs = F,
                       fdr = F,
                       round.to = 3)

map$group.col.spring <- CAF::color.groups(map$groups, cols = c("green",
                                                               "pink",
                                                               "goldenrod"))

map$group.col.fall <- CAF::color.groups(map$groups, cols = c("blue",
                                                             "purple",
                                                             "tomato"))
library(CAF)

res <- tax.shared.multicomp(otu.input = site,
                            var1 = map$season,
                            var2 = map$groups,
                            var1.color1 = map$group.col.spring,
                            var1.color2 = map$group.col.fall,
                            var1.groups = unique( map$season),
                            var2.groups = unique( map$groups),
                            tax.level = tax$family,
                            use.median = F,
                            fdr = F,
                            log = F)
