## coincidence.R

#' Estimate the coincidence function
#'
#' Estimate the coincidence function from backcross data.
#'
#' The coincidence function is the probability of a recombination event in both
#' of two intervals, divided by the product of the two recombination fractions.
#' We estimate this as a function of the distance between the two intervals.
#'
#' Note that we first call \code{\link[qtl]{fill.geno}} to impute any missing
#' genotype data.
#'
#' @param cross Cross object; must be a backcross.  See
#' \code{\link[qtl]{read.cross}} for format details.
#' @param chr Chromosome to consider (only one is allowed).  If missing, the
#' first chromosome is considered.
#' @param pos If provided, these are used as the marker positions.  (This could
#' be useful if you want to do things with respect to physical distance.)
#' @param window Window size used to smooth the estimates.
#' @param fill.method Method used to impute missing data.
#' @param error.prob Genotyping error probability used in imputation of missing
#' data.
#' @param map.function Map function used in imputation of missing data.
#' @return A data.frame containing the distance between intervals and the
#' corresponding estimate of the coincidence.  There are actually two columns
#' of estimates of the coincidence.  In the first estimate, we take a running
#' mean of each of the numerator and denominator and then divide.  In the
#' second estimate, we first take a ratio and then take a running mean.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{gammacoi}}, \code{\link{stahlcoi}}, \code{\link{kfunc}}
#' @references McPeek, M. S. and Speed, T. P. (1995) Modeling interference in
#' genetic recombination.  \emph{Genetics} \bold{139}, 1031--1044.
#' @keywords models
#' @examples
#'
#' map1 <- sim.map(103, n.mar=104, anchor=TRUE, include.x=FALSE, eq=TRUE)
#' x <- sim.cross(map1, n.ind=2000, m=6, type="bc")
#'
#' out <- est.coi(x, window=5)
#' plot(coi1 ~ d, data=out, type="l", lwd=2, col="blue")
#' lines(coi2 ~ d, data=out, lwd=2, col="green")
#' lines(gammacoi(7), lwd=2, col="red", lty=2)
#'
#' @useDynLib xoi
#' @export
est.coi <-
    function(cross, chr, pos, window=0,
             fill.method=c("imp", "argmax"), error.prob=1e-10,
             map.function=c("haldane", "kosambi", "c-f", "morgan"))
{
    if(length(class(cross)) < 2 || class(cross)[1] != "bc")
        stop("This function is only prepared for backcrosses.")

    if(!missing(chr)) {
        if(length(chr) != 1)
            stop("You should specify just one chromosome.")
        cross <- subset(cross, chr=chr)
    }
    else cross <- subset(cross, chr=1)

    dat <- cross$geno[[1]]$data
    if(any(is.na(dat))) {

        fill.method <- match.arg(fill.method)
        map.function <- match.arg(map.function)
        cross <- fill.geno(cross, method=fill.method, error.prob=error.prob,
                           map.function=map.function)
        dat <- cross$geno[[1]]$data
        if(any(is.na(dat))) {
            warning("Some data still missing.")
            dat[is.na(dat)] <- 0
        }
    }

    map <- cross$geno[[1]]$map
    if(!missing(pos)) {
        if(length(pos) != length(map))
            stop("pos must have length ", length(map))
        map <- pos
    }

    ni <- nrow(dat)
    nm <- ncol(dat)
    if(nm < 3) stop("Need at least three markers.")

    npair <- choose(nm-1, 2)

    out <- .C("R_est_coi",
              as.integer(ni),
              as.integer(nm),
              as.integer(npair),
              as.double(map),
              as.integer(dat),
              d=as.double(rep(0, npair)), # distances
              coi1=as.double(rep(0, npair)), # smooth top and bottom and then ratio
              coi2=as.double(rep(0, npair)), # ratio then smooth
              n=as.integer(0), # no. distances to keep
              as.double(window),
              PACKAGE="xoi")

    n <- out$n
    data.frame(d=out$d[1:n], coi1=out$coi1[1:n], coi2=out$coi2[1:n])
}
