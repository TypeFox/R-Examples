library("aroma.core")
library("stats")

x <- matrix(rnorm(1000), ncol=5)
x <- as.data.frame(x)

stats0 <- boxplot(x, plot=FALSE)
stats1 <- lapply(x, FUN=boxplot.stats)
stats1b <- mergeBoxplotStats(stats1)
stopifnot(all.equal(stats0, stats1b))

bxp(stats1b)

