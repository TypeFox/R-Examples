library(tree)
data(cpus, package = "MASS")
cpus$mmin[1:10] <- NA
cpus$mmax[11:nrow(cpus)] <- NA
## segfaulted prior to 1.0-24
try(tree(perf ~ mmin+mmax, cpus))
