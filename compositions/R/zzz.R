#.First.lib <- function(libname,pkgname,...) {
#if( is.null(getOption("robust"))) options(robust=FALSE)
#packageStartupMessage("Welcome to compositions, a package for compositional data analysis.\nFind an intro with \"? compositions\"\n")
#}

.onAttach <- function(libname,pkgname,...) {
if( is.null(getOption("robust"))) options(robust=FALSE)
lapply(objects(pattern="^gsi.*",envir=parent.env(environment())),function(x) if( ! (x %in% objects(pattern="^gsi","package:compositions"))) assign(x,get(x),gsi))
packageStartupMessage("Welcome to compositions, a package for compositional data analysis.\nFind an intro with \"? compositions\"\n")
}

