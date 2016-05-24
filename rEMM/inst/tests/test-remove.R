library("rEMM")
library("testthat")

data("16S")

data <- Mollicutes16S+1
test <- Mollicutes16S[2:10,]+1

## create two EMMs for different data
emm <- EMM("Kullback", threshold=0.1, data=data)
emm


## remove two states
n <- nstates(emm)
e <- remove_clusters(emm, c("3","4"))
n1 <- nstates(e)
expect_equal(n-2L, n1)

## prune w/ copy and w/o copy
rare <- rare_clusters(emm, 5)
e <- prune(emm, 5)
expect_equal(nstates(e), nstates(emm)-length(rare))
expect_equal(clusters(e),setdiff(clusters(emm), rare))


prune(emm, 5, copy=FALSE)
expect_equal(e, emm)

## merge clusters
n <- nclusters(emm)
e <- merge_clusters(emm, c("2", "3", "4", "5"))
expect_equal(nclusters(emm)-3L, nclusters(e))

