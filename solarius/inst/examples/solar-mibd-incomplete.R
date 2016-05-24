
### par
ids <- c("291", "295", "296")

mibddir <- "inst/extdata/solarOutput/solarMibdsCsvIncomplete/"

### data
data(dat30)

dat1 <- dat30

dat2 <- subset(dat1, !(id %in% ids))

ids.fa <- with(dat2, as.character(fa) %in% ids)
ids.mo <- with(dat2, as.character(mo) %in% ids)
ids.par <- ids.fa | ids.mo

dat2 <- within(dat2, {
  fa[ids.par] <- "0"
  mo[ids.par] <- "0"  
})

### models
M1 <- solarMultipoint(trait1 ~ 1, dat1, mibddir = mibddir)
M2 <- solarMultipoint(trait1 ~ 1, dat2, mibddir = mibddir)
