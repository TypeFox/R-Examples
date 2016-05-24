pkgname <- "paran"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('paran')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("paran")
### * paran

flush(stderr()); flush(stdout())

### Name: paran
### Title: Horn's Parallel Analysis of Principal Components/Factors
### Aliases: paran
### Keywords: multivariate

### ** Examples

## perform a standard parallel analysis on the US Arrest data
paran(USArrests, iterations=5000)

## a conservative analysis with different result!
paran(USArrests, iterations=5000, centile=95)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
