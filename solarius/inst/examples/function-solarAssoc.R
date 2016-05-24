data(dat50)

### custom kinship
mod1 <- solarAssoc(traits = "trait", data = phenodata, snpdata = genodata)
mod2 <- solarAssoc(traits = "trait", data = phenodata, snpdata = genodata,
  kinship = kin)

head(sort(mod1$snpf$pSNP))
head(sort(mod2$snpf$pSNP))

### timing when make calculations in parallel
system.time(mod <- solarAssoc(traits = "trait", data = phenodata, snpdata = genodata, cores = 1))

system.time(mod <- solarAssoc(traits = "trait", data = phenodata, snpdata = genodata, cores = 2))


