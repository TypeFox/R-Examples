####  tests/versions.R
## Need this only as long as all other *.R test files, having a matching *.Rout

library("mlmRev", verbose=TRUE)

for(p in c("mlmRev", "lme4", "Matrix")) {
    print(packageDescription(p)); cat(rep.int("-",60), "\n", sep="") }

sessionInfo()
## host specific :
structure(Sys.info()[c(4,5,1:3)], class="simple.list")
