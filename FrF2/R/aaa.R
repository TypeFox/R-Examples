### global variables for avoiding R CMD check notes
### is needed because of file estimable.R
## if (getRversion() >= '2.15.1') globalVariables(".FrF2.currentlychecked")

## create package environment
FrF2Env <- new.env()
getFrF2 <- function(nam) get(nam, envir=FrF2Env) 
putFrF2 <- function(nam, obj) assign(nam, obj, envir=FrF2Env) 