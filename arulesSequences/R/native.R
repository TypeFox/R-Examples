##
## Import native symbols.
##
## ceeboo 2012, 2015

.onLoad <- 
function(libname, pkgname) {
    importNativeSymbol <- function(name, pkg, target = name) {
	obj <- get(name, envir = asNamespace(pkg))
	if (!inherits(obj, "NativeSymbolInfo"))
	    stop(sprintf("'%s' has no native symbol binding in package '%s'",
			 name, pkg)) 
	assign(target, obj, envir = asNamespace(pkgname))
    }

    importNativeSymbol("R_recode_ngCMatrix",    "arules")
    importNativeSymbol("R_asList_ngCMatrix",    "arules")
    importNativeSymbol("R_rowSums_ngCMatrix",   "arules")
    importNativeSymbol("R_colSums_ngCMatrix",   "arules")
    importNativeSymbol("R_cbind_ngCMatrix",     "arules")
    importNativeSymbol("R_rowSubset_ngCMatrix", "arules")
    importNativeSymbol("R_colSubset_ngCMatrix", "arules")
    importNativeSymbol("R_pnindex",             "arules")
    importNativeSymbol("R_pnrindex",            "arules")
}

###
