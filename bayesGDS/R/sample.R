#' @name sample.GDS
#' @aliases sample.GDS
#' @title Collect draws from the target posterior distribution
#' @description Runs the accept-reject phase of the Braun and Damien
#' (2015) algorithm for scalable rejection sampling.
#' @param n.draws number of draws to take from the target posterior
#' density.
#' @param log.phi Vector of log.phi, as computed from the proposal
#' draws.
#' @param post.mode Mode of the target posterior density (numeric
#' vector).
#' @param fn.dens.post Function that returns the log posterior
#' density.  Function should take the parameter vector as the first
#' argument.  Additional arguments are passed as ...
#' @param fn.dens.prop Function that returns the log density of the
#' proposal distribution. The first argument of the function should
#' take either a vector or a matrix.  If the argument is a matrix,
#' each row is considered a sample.  Additional parameters are passed
#' as a list, prop.params.
#' @param fn.draw.prop Function that returns random samples from the
#' proposal density.  This function should return a matrix, with each
#' row being a sample.  Additional parameters are passed as a list,
#' prop.params.
#' @param prop.params Object (list or vector) to be passed to both
#' fn.dens.prop and fn.draw.prop.Contains parameters for the proposal
#' distribution.  See details.
#' @param ... Additional parameters to be passed to fn.dens.post.
#' @param max.tries Maximum number of proposal draws to try, without a
#' success.  This prevents the routine from being stuck in an endless
#' loop.
#' @param report.freq The frequency that the function will report the
#' current iteration.  For example, if report.freq=5, the function
#' will display a message after every fifth iteration.
#' @param announce If TRUE, will print a message when a proposal is
#' accepted as a sample from the target posterior distribution.
#' @param thread.id An identifier used in the announce function.  This
#' is useful if running sample.GDS on multiple processors, to collect
#' multiple batches of samples. Defaults to 1.
#' @param seed Sets a random seed within the call to sample.GDS.
#' Useful for assigning different seeds to calls to sample.GDS that
#' are running on different threads or processors.  Defaults to
#' .Random.seed.
#' @return a list with the following elements:
#' \item{draws}{A matrix with each draw in a row, and each parameter
#' in a column}
#' \item{counts}{The number of attempts that it took to get an
#' accepted draw.  The accepted draw counts, so the count will always
#' be at least 1.}
#' \item{gt.1}{A vector that indicates if the phi for that draw was
#' greater than 1.  Available as a diagnostic.  Normally, these should
#' all be FALSE.  Any values of TRUE suggest that a change in proposal
#' density might be warranted.}
#' \item{log.post.dens}{A numeric vector.  Log posterior density for
#' each draw.}
#' \item{log.prop.dens}{A numeric vector. Log of the proposal density
#' for each draw.}
#' \item{log.thresholds}{Vector of threshold draws (log u) from the
#' accept-reject algorithm.  Sorted in ascending order.}
#' \item{log.phi}{A numeric vector.  Value of log.phi for the accepted
#' draws.}
#' @references
#' Braun, Michael and Paul Damien (2015).  Scalable Rejection Sampling for
#' Bayesian Hierarchical Models. Marketing Science. Articles in Advance.
#' http://doi.org/10.1287/mksc.2014.0901
#' @export
sample.GDS <- function(n.draws, log.phi, post.mode, fn.dens.post,
                       fn.dens.prop, fn.draw.prop,
                       prop.params, ...,
                       max.tries=1000000, report.freq=1,
                       announce=FALSE, thread.id=1, seed=.Random.seed) {
    if (any(!is.finite(seed))) {
        stop("Error in sample.GDS:  all values in seed must be finite")
    }
    set.seed(seed)
    v <- get.cutoffs(log.phi, n.draws)
    nvars <- length(post.mode)
    log.c1 <- fn.dens.post(post.mode, ...)
    log.c2 <- fn.dens.prop(post.mode, prop.params)

    v.keep <- vector("numeric",length=n.draws)
    counts <- vector("integer",length=n.draws)
    gt.1 <- vector("integer",length=n.draws)
    draws.keep <- matrix(nrow=n.draws, ncol=nvars)
    log.post.dens <- vector("numeric",length=n.draws)
    log.prop.dens <- vector("numeric",length=n.draws)
    crit.keep <- vector("numeric",length=n.draws)

    remaining.draws <- n.draws
    idx.1 <- 1
    count.idx <- 1

    while((remaining.draws>0) & (count.idx <= max.tries)) {

        if ((count.idx %% report.freq) == 0) {
            cat("thread ",thread.id,"  count ",count.idx,
                "  remaining draws = ",remaining.draws,"\n")
        }

        x <- as.matrix(fn.draw.prop(remaining.draws, prop.params))
        dens.1 <- apply(x,1,fn.dens.post,...)
        dens.1 <- as.vector(dens.1)
        dens.2 <- fn.dens.prop(x,prop.params)
        dens.2 <- as.vector(dens.2)
        crit <- dens.1-dens.2-log.c1+log.c2
        ww <- which(-v < crit)
        n.keep <- length(ww)
        if (any(crit>0)) {
            cat("bayesGDS warning:  accepted draw(s) with log.phi>0\n")
        }

        if (n.keep > 0) {
            idx.range <- idx.1:(idx.1 + n.keep - 1)
            v.keep[idx.range] <- v[ww]
            counts[idx.range] <- count.idx
            draws.keep[idx.range,] <- x[ww,]
            log.post.dens[idx.range] <- dens.1[ww]
            log.prop.dens[idx.range] <- dens.2[ww]
            gt.1[idx.range] <- crit[ww]>0
            crit.keep[idx.range] <- crit[ww]
            idx.1 <- idx.1+n.keep
            remaining.draws <- remaining.draws-n.keep
            v <- v[-ww] ## drop v that were successful
            if (announce) {
                cat("Sample accepted in thread ",thread.id,". count =  ",count.idx,"\n")
            }
        }
        count.idx <- count.idx+1
    }

    res <- list(draws=draws.keep,
                counts=counts,
                gt.1=gt.1,
                log.post.dens=log.post.dens,
                log.prop.dens=log.prop.dens,
                log.thresholds=v.keep,
                log.phi=crit.keep)

    return(res)
}
