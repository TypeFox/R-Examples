# Test the C implementation of ranking entropy (about 40x speed up)
library(smaa)

data <- dget('rankEntropy.txt')
N <- dim(data$ranks)[1]
m <- dim(data$ranks)[2]

counts <- .Call("smaa_countRankings", t(data$ranks))

stopifnot(all.equal(counts[counts > 0], data$counts))
stopifnot(all.equal(smaa.entropy.ranking(data$ranks), data$entropy))
