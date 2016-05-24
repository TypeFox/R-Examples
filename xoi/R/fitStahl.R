## fitStahl.R

#' Calculate log likelihood for Stahl model
#'
#' Calculate the log likelihood for the Stahl model for varying parameters,
#' with data on crossover locations.
#'
#' See Housworth and Stahl (2003) and Broman and Weber (2000) for details of
#' the method.
#'
#' If neither \code{nu} nor \code{p} has length 1, they both must have the same
#' length.  If one has length 1 and the other does not, the one with length 1
#' is repeated so that they both have the same length.
#'
#' @param xoloc A list of crossover locations (in cM), each component being a
#' vector of locations for a different meiotic product.
#' @param chrlen Chromosome length (in cM), either of length 1 or the same
#' length as \code{xoloc}.
#' @param nu A vector of interference parameters (\eqn{\nu}{nu}) at which to
#' calculate the log likelihood.
#' @param p A vector of parameter values for the proportion of crossovers from
#' the no interference pathway.
#' @param max.conv Maximum limit for summation in the convolutions to get
#' inter-crossover distance distribution from the inter-chiasma distance
#' distributions.  This should be greater than the maximum number of chiasmata
#' on the 4-strand bundle.
#' @param integr.tol Tolerance for convergence of numerical integration.
#' @param max.subd Maximum number of subdivisions in numerical integration.
#' @param min.subd Minimum number of subdivisions in numerical integration.
#' @return A vector of log likelihoods.
#'
#' The corresponding values of nu and p are saved as attributes.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link[qtl]{fitstahl}}
#' @references Housworth, E. A. and Stahl, F. W. (2003) Crossover interference
#' in humans. \emph{Am. J. Hum. Genet.} \bold{73}, 188--197.
#'
#' Broman, K. W. and Weber, J. L. (2000) Characterization of human crossover
#' interference. \emph{Am. J. Hum. Genet.} \bold{66}, 1911--1926.
#' @keywords models
#' @examples
#'
#' data(bssbsb)
#' xoloc <- find.breaks(bssbsb, chr=1)
#'
#' loglik <- stahlLoglik(xoloc, nu=4, p=c(0.05, 0.1, 0.15))
#'
#' @useDynLib xoi
#' @export
stahlLoglik <-
    function(xoloc, chrlen, nu, p,
             max.conv=25, integr.tol=1e-8, max.subd=1000, min.subd=10)
{
    if(is.data.frame(xoloc)) stop("xoloc should not be a data.frame.")
    if(!is.list(xoloc)) stop("xoloc should be a list.")

    if(missing(chrlen) && "L" %in% names(attributes(xoloc)))
        chrlen <- attr(xoloc, "L")

    if(length(chrlen) == 1) {
        chrlen <- rep(chrlen, length(xoloc))
        constant.chrlen <- TRUE
    }
    else {
        constant.chrlen <- FALSE
        if(length(chrlen) != length(xoloc))
            stop("chrlen should have length 1 or the same as length(xoloc).")
    }

    if(length(nu)==1) {
        if(length(p) > 1) nu <- rep(nu, length(p))
    }
    else {
        if(length(p)==1) p <- rep(p, length(nu))
        else if(length(p) != length(nu))
            stop("nu and p should be the same length (though either can have length 1).")
    }

    flag <- 0
    for(i in seq(along=xoloc)) {
        thisxoloc <- unlist(xoloc[[i]])
        if(any(thisxoloc < 0 | thisxoloc > chrlen[i])) {
            flag <- 1
            break
        }
    }
    if(flag) stop("xoloc should be between 0 and chrlen.")

    # intercross?  then send to stahlLoglikF2
    if(is.list(xoloc[[1]]))
        return(stahlLoglikF2(xoloc, chrlen, nu, p, max.conv, integr.tol, max.subd, min.subd, constant.chrlen))

    n <- sapply(xoloc, length)
    n.nu <- length(nu)

    loglik <- .C("R_stahl_loglik",
                 as.integer(length(xoloc)),
                 as.integer(n),
                 as.double(unlist(xoloc)/100), # convert to Morgans
                 as.double(chrlen/100), # convert to Morgasn
                 as.integer(n.nu),
                 as.double(nu),
                 as.double(p),
                 loglik=as.double(rep(0,n.nu)),
                 as.integer(max.conv),
                 as.double(integr.tol),
                 as.integer(max.subd),
                 as.integer(min.subd),
                 as.integer(constant.chrlen),
                 PACKAGE="xoi")$loglik
    attr(loglik, "nu") <- nu
    attr(loglik, "p") <- p

    loglik
}

# stahlLoglik for intercross
stahlLoglikF2 <-
    function(xoloc, chrlen, nu, p,
             max.conv=25, integr.tol=1e-8, max.subd=1000, min.subd=10,
             constant.chrlen=FALSE)
{
    n.ind <- length(xoloc)
    n.alternatives <- unlist(lapply(xoloc, length))
    n.xo.per <- unlist(lapply(xoloc, lapply, lapply, length))
    n.products <- length(n.xo.per) # = 2 * sum(n.alternatives)

    # expand chromosome lengths
    chrlen <- rep(chrlen, n.alternatives*2)

    n.nu <- length(nu)

    loglik <- .C("R_stahl_loglik_F2",
                 as.integer(n.ind),
                 as.integer(n.alternatives),
                 as.integer(n.products),
                 as.integer(n.xo.per),
                 as.double(unlist(xoloc)/100), # convert to Morgans
                 as.double(chrlen/100), # convert to Morgans
                 as.integer(n.nu),
                 as.double(nu),
                 as.double(p),
                 loglik=as.double(rep(0,n.nu)),
                 as.integer(max.conv),
                 as.double(integr.tol),
                 as.integer(max.subd),
                 as.integer(min.subd),
                 as.integer(constant.chrlen),
                 PACKAGE="xoi")$loglik
    attr(loglik, "nu") <- nu
    attr(loglik, "p") <- p

    loglik
}

######################################################################
# the same, but with the arguments reordered and with p and nu stuck
# together, and assuming they have length 1
######################################################################
fitStahl.sub <-
    function(param, xoloc, chrlen, max.conv=25, integr.tol=1e-8,
             max.subd=1000, min.subd=10)
{
    if(param[1] < 0 || param[2] < 0 || param[2] > 1)
        return(Inf)

    -stahlLoglik(xoloc, chrlen, param[1], param[2], max.conv,
                 integr.tol, max.subd, min.subd)
}

# here, to optimize for the model with p=0
fitStahl.sub2 <-
    function(nu, xoloc, chrlen, max.conv=25, integr.tol=1e-8,
             max.subd=1000, min.subd=10)
{
    if(nu < 0) return(Inf)

    -stahlLoglik(xoloc, chrlen, nu, 0, max.conv,
                 integr.tol, max.subd, min.subd)
}

# function to optimize for the Stahl model


#' Fit Stahl model
#'
#' Fit the Stahl model for crossover interference to data on crossover
#' locations.
#'
#' See Housworth and Stahl (2003) and Broman and Weber (2000) for details of
#' the method.
#'
#' We first use \code{\link[stats]{optimize}} to find the MLE with the
#' contraint \code{p=0}, followed by use of \code{\link[stats]{optim}} to do a
#' 2-dimensional optimization for the MLEs of the pair.
#'
#' @param xoloc A list of crossover locations (in cM), each component being a
#' vector of locations for a different meiotic product.
#' @param chrlen Chromosome length (in cM), either of length 1 or the same
#' length as \code{xoloc}.
#' @param nu Interference parameter (\eqn{\nu}{nu}).  This should be a pair of
#' values to be used as endpoints to first do a 1-dimensional optimization with
#' \eqn{p=0}.
#' @param p Starting value for the proportion of crossovers from the no
#' interference pathway, for the 2-dimensional optimization.
#' @param max.conv Maximum limit for summation in the convolutions to get
#' inter-crossover distance distribution from the inter-chiasma distance
#' distributions.  This should be greater than the maximum number of chiasmata
#' on the 4-strand bundle.
#' @param integr.tol Tolerance for convergence of numerical integration.
#' @param max.subd Maximum number of subdivisions in numerical integration.
#' @param min.subd Minimum number of subdivisions in numerical integration.
#' @param verbose If TRUE, print tracing information.  If "\dots{}" includes
#' \code{control}, this is ignored.
#' @param \dots Further arguments sent to \code{\link[stats]{optim}}.
#' @return A vector with the estimates of \eqn{\nu}{nu} (interference
#' parameter) and \eqn{p} (proportion of crossovers coming from the no
#' interference pathway), the maximized log likelihood, the estimate of nu with
#' p constrained to be 0, the maximized log likelihood in this case, and the
#' log likelihood ratio for comparing the model with p allowed to vary freely
#' versus contrained to be 0.  (Note that it's the natural log of the
#' likelihood ratio, and not twice that.)
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{fitGamma}}, \code{\link{stahlLoglik}},
#' \code{\link{simStahl}}
#' @references Housworth, E. A. and Stahl, F. W. (2003) Crossover interference
#' in humans. \emph{Am. J. Hum. Genet.} \bold{73}, 188--197.
#'
#' Broman, K. W. and Weber, J. L. (2000) Characterization of human crossover
#' interference. \emph{Am. J. Hum. Genet.} \bold{66}, 1911--1926.
#' @keywords models
#' @examples
#'
#' data(bssbsb)
#' \dontshow{bssbsb <- bssbsb[,1:50]}
#'
#' xoloc <- find.breaks(bssbsb, chr=1)
#' L <- attr(xoloc, "L")
#'
#' # get MLE (limiting maximum iterations to 10, just for speed in this example)
#' \dontrun{mle <- fitStahl(xoloc, L, nu=c(9, 12), control=list(maxit=10))}
#' \dontshow{mle <- fitStahl(xoloc, L, nu=c(9, 12), control=list(maxit=2))}
#'
#' @importFrom stats optimize optim
#' @export
fitStahl <-
    function(xoloc, chrlen, nu=c(1,20), p=0.02, max.conv=25, integr.tol=1e-8,
             max.subd=1000, min.subd=10, verbose=TRUE, ...)
{
    if(is.data.frame(xoloc)) stop("xoloc should not be a data.frame.")
    if(!is.list(xoloc)) stop("xoloc should be a list.")

    if(missing(chrlen) && "L" %in% names(attributes(xoloc)))
        chrlen <- attr(xoloc, "L")

    if(length(nu) > 2) {
        warning("nu should have length 2; using the first two values.")
        nu <- nu[1:2]
    }
    if(length(nu) != 2)
        stop("nu should have length 2.")
    if(length(p) > 1) {
        warning("p should have length 1; using the first value.")
        p <- p[1]
    }
    if(length(p) != 1)
        stop("p should have length 1.")

    out0 <- optimize(fitStahl.sub2, interval=c(nu[1], nu[2]), xoloc=xoloc, chrlen=chrlen,
                     max.conv=max.conv, integr.tol=integr.tol, max.subd=max.subd,
                     min.subd=min.subd)
    nu <- out0$minimum
    if(verbose) cat("For p=0, nuhat =", nu, "\n       log lik =", -out0$objective, "\n")

    if(verbose>1 && !("control" %in% names(list(...))))
        out <- optim(c(nu, p), fitStahl.sub, xoloc=xoloc, chrlen=chrlen,
                     max.conv=max.conv, integr.tol=integr.tol, max.subd=max.subd,
                     min.subd=min.subd, control=list(trace=verbose-1), ...)
    else
        out <- optim(c(nu, p), fitStahl.sub, xoloc=xoloc, chrlen=chrlen,
                     max.conv=max.conv, integr.tol=integr.tol, max.subd=max.subd,
                     min.subd=min.subd, ...)

    if(verbose)
        cat("\n  nuhat =", out$par[1], "\n   phat =", out$par[2], "\nlog lik =", -out$value, "\n")

    if(out0$objective <= out$value) {
        out <- c(out0$minimum, 0, -out0$objective, out0$minimum, -out0$objective, 0)
        if(verbose) cat("Inferred that p=0\n")
    }
    else {
        out <- c(out$par, -out$value, out0$minimum, -out0$objective, out0$objective - out$value)
        if(verbose) cat("Inferred that p>0\n")
    }
    names(out) <- c("nu", "p", "loglik", "nu0", "loglik0", "ln LR testing p=0")
    out
}
