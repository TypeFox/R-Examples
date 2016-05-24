##
## Generic for extractor of x-variables
setGeneric("cps", function(cpd, ctrl) standardGeneric("cps"))
##
## Generic for extractor of x-variables
setGeneric("getx", function(object) standardGeneric("getx"))
##
## Generic for extractor of y-variables
setGeneric("gety", function(object) standardGeneric("gety"))
##
## Generic for extractor of s-variables
setGeneric("gets", function(object) standardGeneric("gets"))
##
## Generic for extractor of z-variables
setGeneric("getz", function(object) standardGeneric("getz"))
##
## Generic for extractor of state of convex program
setGeneric("getstate", function(object) standardGeneric("getstate"))
##
## Generic for extractor of optimizer's status
setGeneric("getstatus", function(object) standardGeneric("getstatus"))
##
## Generic for extractor of number of iterations
setGeneric("getniter", function(object) standardGeneric("getniter"))
##
## Generic for extractor of control parameters
setGeneric("getparams", function(object) standardGeneric("getparams"))
