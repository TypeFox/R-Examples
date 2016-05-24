# This function returns the solution of the minimization problem
setGeneric("zsummin",
           function (o, x=0, y=0, lp=numeric(0), max.iter=100, eps=1.e-3, verbose=FALSE, algorithm="weiszfeld", ...) standardGeneric("zsummin")
)

# Check lp value and call to the specific function
setMethod("zsummin", "loca.p",
function (o, x=0, y=0, lp=numeric(0), max.iter=100, eps=1.e-3, verbose=FALSE, algorithm="weiszfeld", ...)
   {
     if (length(lp) == 0) return(zsuml2min(o=o, x=x, y=y, max.iter=max.iter, eps=eps, verbose=verbose, algorithm=algorithm, ...))
     else if (lp >= 1) return(zsumlpmin(o=o, x=x, y=y, p=lp, max.iter=max.iter, eps=eps, verbose=verbose, algorithm=algorithm, ...))
     else stop(paste(lp, gettext("is not a valid value for lp, use 1 <= lp", domain = "R-orloca")))
   }
)

warning.max.iter <- function(max.iter)
  {
    warning(paste(gettext("zsummin: Maximun number of iteration reached", domain = "R-orloca"), " (max.iter = ", max.iter, ")\n", gettext("The solution may be non-optimal.\nPerhaps, you could try increasing max.iter value.", domain = "R-orloca"), sep=""), call. = F)
  }
