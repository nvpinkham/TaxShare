
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
