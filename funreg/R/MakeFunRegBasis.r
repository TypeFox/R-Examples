#'@title Make basis for functional regression (for internal use by other package functions)
#' @description This is a function for internal use (i.e., a user will not need 
#' to call it directly for usual data analysis tasks).  Recall that
#' functional coefficients are estimated as a linear combination of 
#' basis functions, thus changing a nonparametric into a parametric
#' estimation problem.  This function constructs the matrix of basis 
#' function values for doing a functional regression. 
#' @param basis.type is a character string, either \code{TruncatedPower}
#' or \code{BSpline}. This tells whether the basis functions should be
#' calculated as B-splines (see Eilers and Marx, 1996) or 
#' as truncated power splines (see Ruppert, Wand, and Carroll, 2003).
#' @param deg is the degree of the basis functions (roughly, their amount 
#' of complexity) and should generally be 1, 2, or 3. 
#' @param num.knots is the number of knots in the basis; the higher 
#' this is, the more flexible the estimated function will be.  
#' If it is too low, the estimated function may be too simple
#'  (i.e.,biased towards being too smooth).  If it is
#' too high, the function may be hard to interpret.
#' @param times is the vector of measurement times (more technically, 
#' real-valued index values for the functional covariate) at which
#' the basis functions should be evaluated. 
#' @return Returns a list with two components.  The first,
#' \code{interior.knot.locations}, tells the selected locations on the time
#' axis for each interior knot. The second, \code{basis.for.betafn},
#' is a matrix with one row for each time value in the input vector
#' \code{times} and one column for each basis function. It represents
#' the values of the basis functions themselves.   
#' @references Eilers, P. H. C., and Marx, B. D. (1996). 
#' Flexible smoothing with B-splines and penalties (with
#' comments and rejoinder). Statistical Science, 11, 89-121.
#' 
#'  Ruppert, D., Wand, M. P., and Carroll,
#'  R. J. (2003) Semiparametric regression. Cambridge: Cambridge.
#'@importFrom splines spline.des
#'@export
make.funreg.basis <- function( basis.type,
                                deg,
                                num.knots,
                                times ) {
    # An internal helper function to make the basis for the
    # spline approximation of the beta function.
    if (basis.type=="TruncatedPower") {
        # Follows Ruppert, Wand and Carroll (2003)
        num.intervals <- num.knots + 1;
        dx <- (max(times)-min(times))/num.intervals;
        all.knot.locations <- seq(min(times)+dx,
                                  max(times)-dx,
                                  by=dx);
        interior.knot.locations <- all.knot.locations;
        basis.for.betafn <- NULL;
        for (d in 0:deg) {
            basis.for.betafn <- cbind(basis.for.betafn,
                                    times^d);
        }
        basis.for.betafn <- cbind(basis.for.betafn,
                                pmax(outer(times,all.knot.locations,"-"),0)^deg);
    }
    else if (basis.type=="BSpline") {
        # Follows Eilers & Marx (1996)
        num.intervals <- num.knots-2*deg;
        dx <- (max(times)-min(times))/num.intervals;
        all.knot.locations <- seq(min(times)-deg*dx,
                                  max(times)+deg*dx,
                                  by=dx);
        interior.knot.locations <- seq(min(times),
                                       max(times),
                                       by=dx);
        basis.for.betafn <- spline.des(all.knot.locations,
                                     times,
                                     deg+1,
                                     rep(0,length(times)),
                                     outer.ok=TRUE)$design;
    }
    else stop();
    invisible(list(interior.knot.locations=interior.knot.locations,
                   basis.for.betafn=basis.for.betafn));
}