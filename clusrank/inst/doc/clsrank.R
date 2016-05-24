## ----rank sum test for unbalanced data-----------------------------------
library(clusrank)
data(crd)
## Using large-sample test.
cluswilcox.test(z ~ group(group) + cluster(id), data = crd)
## For small sample, using the permutation test.
cluswilcox.test(z ~ group(group) + cluster(id), data = crd, permutation = TRUE)

## ----rank sum test for stratified data-----------------------------------
data(crdStr)
cluswilcox.test(z ~ group(group) + cluster(id) + stratum(stratum), data = crdStr)

## ----signed rank test for balanced data----------------------------------
## Large sample signed rank test for clustered data
data(crsd)
cluswilcox.test(z, cluster = id, data = crsd)
## Small sample test
data(crsd)
cluswilcox.test(z, cluster = id, data = crsd, permutation = TRUE)

## ----signed rank test for unbalanced data--------------------------------
data(crsdUnb)
cluswilcox.test(z, cluster = id, data = crsdUnb)

## ----CARMS1234-----------------------------------------------------------
## Carry out clustered rank sum test for the subset
## with CARMS grade 1, 2, 3 and 4.
data(sedlab)
cluswilcox.test(CARMS ~ cluster(ID) + stratum(Agesex) + group(Variant),
                data = sedlab, subset = CARMS %in% c(1, 2, 3, 4))

## ----CARMS1235-----------------------------------------------------------
## Carry out clustered rank sum test for the subset
## with CARMS grade 1, 2, 3 and 5.
data(sedlab)
cluswilcox.test(CARMS ~ cluster(ID) + stratum(Agesex) + group(Variant),
                data = sedlab, subset = CARMS %in% c(1, 2, 3, 5))

