require("pdc")

# convert warnings to errors
options(warn=2)

pdcDist(X=matrix(rnorm(100),ncol=2) )
