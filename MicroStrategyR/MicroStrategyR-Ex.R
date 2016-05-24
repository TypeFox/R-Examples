pkgname <- "MicroStrategyR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('MicroStrategyR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("deployR")
### * deployR

flush(stderr()); flush(stdout())

### Name: deployR
### Title: deployR Function
### Aliases: deployR
### Keywords: microstrategy

### ** Examples

## Not run: 
##D deployR()
## End(Not run)



### * <FOOTER>
###
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
