
source("R/Functions.2026.01.18.R")

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
                              use.median = T)
par(mar = c(5, 10, 4, 4))
res <- tax.shared.plot(tax.res = tax.table,
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
                               use.median = F,
                               round.to = 3)

res2 <- tax.shared.plot(tax.res = tax.table2,
                        color1 = "salmon1",
                        color2 = "lightgreen",
                        log = F,
                        round.to = 3)

map$group.col.spring <- CAF::color.groups(map$groups, cols = c("green",
                                                               "pink",
                                                               "goldenrod"))

map$group.col.fall <- CAF::color.groups(map$groups, cols = c("blue",
                                                               "purple",
                                                               "tomato"))

res <- tax.shared.multicomp(otu.input = site,
                            var1 = map$season,
                            var2 = map$groups,
                            var1.color1 = map$group.col.spring,
                            var1.color2 = map$group.col.fall,
                            var1.groups = unique( map$season),
                            var2.groups = unique( map$groups),
                            tax.level = tax$family,
                            use.median = T,
                            family.name.length = 3,
                            fdr = F,
                            log = F)

map <- read.csv("data/map.csv", row.names = 1)
otu <- read.csv("data/otu.csv", row.names = 1)
tax <- read.csv("data/tax.csv")

fam <- aggregate(t(otu), list(tax$family), sum)
fam <- fam[order(rowSums(fam[,-1]), decreasing = T) , ]

p <- map$trial == "Huston ABX Trial 1" &  map$antibiotic.group == "Ibezapolstat"
p <- map$antibiotic.group == "Ibezapolstat"

otu <- as.data.frame(otu)

#otu <- otu[ , tax$family == "Odoribacteraceae"|tax$family == "Akkermansiaceae"]
#tax <- tax[tax$family == "Odoribacteraceae"|tax$family == "Akkermansiaceae" , ]

tax.anti <- tax.shared.table(otu.input = otu[p,],
                               var =map$day2anti[p],
                               group1 = "0",
                               group2 = "10",
                               min.obs = 1,
                               tax.level = tax$family,
                               PresAbs = F,
                               fdr = F,
                               use.median = F,
                               round.to = 3)

par(mar = c(5, 30, 4, 4))

a <- colSums(tax.anti)
as.numeric(a) > 0

tax.anti.pick <- tax.anti[, colSums(tax.anti) > 0]

res <- tax.shared.plot(tax.res = tax.anti.pick,
                       log = T, family.name.length = 3)

p <- map$day2anti %in% c(0, 10) & map$trial == "Huston ABX Trial 1"
p <- map$day2anti %in% c(0, 10)

fams <- aggregate(t(otu[p,]), list(tax$family), median)
fams <- fams[order(rowSums(fams[,-1]), decreasing = T) , ]

unique(var1.color2)

tax.shared.multicomp(otu.input = otu [p,],
                     var1 = map$day2anti [p],
                     var2 = map$antibiotic.group [p],
                     var1.groups =c("0", "10") ,
                     var2.groups = unique(map$antibiotic.group),
                     tax.level = tax$family ,
                     var1.color1 = "grey" ,
                     var1.color2 = map$col.group ,
                     PresAbs = F ,
                     log = T ,
                     use.median = F,
                     fdr = F )

p <- map$day2anti %in% c(0, 10) & map$antibiotic.group != "ND Control"
par(mar = c(5, 20, 4, 20))

tax.shared.multicomp(otu.input = otu [p,],
                     var1 = map$day2anti [p],
                     var2 = map$antibiotic.group [p],
                     var1.groups =c("0", "10") ,
                     var2.groups = unique(map$antibiotic.group[p]),
                     tax.level = tax$family ,
                     var1.color1 = "grey" ,
                     var1.color2 = map$col.group[p] ,
                     PresAbs = F ,
                     log = T ,
                     use.median = F,
                     family.name.length = 3,
                     taxa2include = fams$Group.1[1:5],
                     fdr = F )



