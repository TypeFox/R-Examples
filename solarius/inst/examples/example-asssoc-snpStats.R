library(snpStats)
data(families)

#tests <- tdt.snp(data = pedData, snp.data = genotypes)

phen <- within(pedData, {
  id <- rownames(pedData)
})
phen <- rename(phen, c(familyid = "famid"))

mod <- solarPolygenic(affected ~ 1, phen, dir = "solar")
