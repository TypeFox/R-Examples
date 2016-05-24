pkgname <- "regtest"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('regtest')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("binregtest")
### * binregtest

flush(stderr()); flush(stdout())

### Name: binregtest
### Title: Binary regression test
### Aliases: binregtest
### Keywords: debugging documentation

### ** Examples

wronglog <- function(x, base=exp(1)){
  if (x>0)
    log(x, base=base)
  else
    NA
}
binregtest(wronglog, log, x=as.list(0:3), base=list(2, exp(1), 10))



cleanEx()
nameEx("is.all.equal")
### * is.all.equal

flush(stderr()); flush(stdout())

### Name: is.all.equal
### Title: wrapper for all.equal
### Aliases: is.all.equal
### Keywords: debugging utilities

### ** Examples

  all.equal(1,2)
  is.all.equal(1,2)



cleanEx()
nameEx("timefactor")
### * timefactor

flush(stderr()); flush(stdout())

### Name: timefactor
### Title: compare timing of two expressions
### Aliases: timefactor
### Keywords: misc

### ** Examples

  timefactor(Sys.sleep(0.1), Sys.sleep(1), 10, 1)



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
