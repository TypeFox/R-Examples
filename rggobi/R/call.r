
# GGobi symbol creation
# Maps the given name to the name of the corresponding C routine
#
# A simple way of generating the prefix for a symbol
# used in this package/library so that we can hide
# it from other packages and avoid conflicts.
#
# @keyword dynamic 
# @keyword internal
# @value the name of the C routine corresponding to its argument
.ggobi.symbol <- function(name) paste("RS_GGOBI", name, sep="_")

# Calling native routines
# Wrappers for calling C routines in the R-ggobi library.
# 
# \code{.GGobiC} and \code{.GGobiCall} convert the name and then call 
# their C invocation counterparts.
# 
# These functions map the simple name of a C routine into the
# package-specific version of that name.  These allow use to hide the
# use a name \emph{mangling} scheme of our choosing for the C level
# routines in the shared library/DLL that provides the glue between R
# and ggobi.  This is useful for avoiding name conflicts with other C
# code in R or other packages.  These are only of relevance to the
# developers of this package and those working with its C code.
# 
# The mapping of the name to its corresponding C routine name
# is done in conjunction with the pre-processor macro
# \code{RS_GGOBI}. These  must be synchronized.
# 
# @alias .GGobiC
# @arguments  the simple name of the C routine to be resolved
# @arguments the arguments that to be passed to the \code{\link{.C}} or \code{\link{.Call}}
# @arguments the ggobi instance identifier that is to be passed to the C routine as its last argument
# @value the same result as the corresponding \code{.C} and \code{.Call}
# @references \url{http://www.ggobi.org/}
# @seealso \code{\link{.C}}, \code{\link{.Call}}
# @keyword dynamic 
# @keyword internal
.GGobiCall <- function(.name, ..., .gobi = ggobi_get(), .test=TRUE) {
	if (.test && !is.null(.gobi) && !valid_ggobi(.gobi)) stop("Invalid ggobi reference", call.=FALSE)
	.Call(.ggobi.symbol(.name), ..., .gobi, PACKAGE = "rggobi")
}

.GGobiC <- function(.name, ..., .gobi = ggobi_get(), .test=TRUE) {
	if (.test && !is.null(.gobi) && !valid_ggobi(.gobi)) stop("Invalid ggobi reference", call.=FALSE)
        sym <- .ggobi.symbol(.name)
        if (!is.null(.gobi))
          .C(sym, ..., .gobi, PACKAGE = "rggobi")
        else .C(sym, ..., PACKAGE = "rggobi")
}

# Validity checking
# Determines whether a reference to an internal ggobi object is valid
#
# One can create multiple, independent ggobi instances within a single
# R session and one can also remove them either programmatically or
# via the GUI.  To be able to refer to these objects which are
# actually C-level internal objects, one has a reference or handle
# from an S object. Since the C level object can be destroyed while the S 
# object still refers to them, this function allows one to check whether the 
# internal object to which R refers is still in existence.
#
# @arguments an object of class \code{ggobi} which refers to an internal ggobi instance.
# @value \code{TRUE} if real object still exist, \code{FALSE} otherwise
# @keyword dynamic 
#X g <- ggobi(mtcars)
#X valid_ggobi(g)
#X close(g)
#X valid_ggobi(g) 
# @keyword internal
valid_ggobi <- function(.gobi) {
	.GGobiCall("isValid", .gobi=.gobi, .test=FALSE)
}
