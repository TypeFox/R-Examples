# This tests part of a problem as of zoo 1.7-11 that is fixed in the devel version.
#changeTSrepresentation() seems sometimes to need  require("zoo") rather than
# requireNamespace("zoo"). zoo needs to overwrite the as.Date() generic since there 
# is no  default for the origin in the as.Date.numeric() method.

#test2 <- function(){
#     loadNamespace("zoo")
#     z0 <- ts(matrix(rnorm(10),10,1), start=c(1990,1), frequency=1)
#     z1 <- zoo::as.Date(stats::time(z0))
#     as.ts(z1)
#     }

require("tframePlus")

z0 <- ts(matrix(rnorm(10),10,1), start=c(1990,1), frequency=1)

z <- changeTSrepresentation(z0, "zoo")

z <- changeTSrepresentation(z0, "ts")

z

if(start(z) != 1990 ) stop("zoo changeTSrepresentation test 1 failed.")

if(frequency(z) != 1) stop("zoo changeTSrepresentation test 2 failed.")

if( ! z == z0)        stop("zoo changeTSrepresentation test 3 failed.")

if(class(z) !=  "ts") stop("zoo changeTSrepresentation test4 failed.")
