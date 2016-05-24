## kwak_coincidence.R

#' Estimate coincidence function
#'
#' Estimate coincidence function for a chromosome.
#'
#' @param cross Cross object; must be a backcross.  See
#' \code{\link[qtl]{read.cross}} for format details.
#' @param chr Chromosome to consider (only one is allowed).  If missing, the
#' first chromosome is considered.
#' @param window Window size
#' @param ncalc Total number of points for calculations.
#' @return Data frame with columns \code{distance} and \code{coincidence}.  The
#' input argument \code{window} is kept as an attribute.
#' @author Il youp Kwak
#' @seealso \code{\link{intensity}}, \code{\link{est.coi}}
#' @keywords utilities
#' @examples
#'
#' map1 <- sim.map(103, n.mar=104, anchor=TRUE, include.x=FALSE, eq=TRUE)
#' x <- sim.cross(map1, n.ind=2000, m=6, type="bc")
#'
#' out <- coincidence(x, ncalc=101)
#' plot(out, type="l", lwd=2, ylim=c(0, max(out[,2])))
#'
#' @useDynLib xoi
#' @export
coincidence <-
    function(cross, chr, window=5, ncalc=500)
{
    if(!missing(chr)) {
        cross <- subset(cross, chr)
        if(nchr(cross) > 1)
            warning("Considering only chr ", names(cross$geno)[1])
    }

    if(class(cross)[1] != "bc")
        stop("coincidence() currently working only for a backcross.")

    g <- cross$geno[[1]]$data
    g[is.na(g)] <- 0
    map <- cross$geno[[1]]$map

    xomat <- identify_xo(g)$xomat
    n.ind <- nind(cross)

    coincidenceSUB(xomat, window, map, n.ind, ncalc)
}

coincidenceSUB <-
    function(xomat, window, marker, n_ind, N)
{

    inten <- intensitySUB(xomat, window, marker, n_ind, N)
    int_dat <- inten$intensity*window/100
    center = inten$position

    n_pos <- length(marker)
    n_center <- length(center)
    n_xo <- ncol(xomat)

    start_d <- min(which(center >= window/2))-1

    output <- .C("R_get_coincidence",
                 as.integer(xomat),
                 as.double(int_dat),
                 as.double(window),
                 as.double(center),
                 as.integer(n_xo),
                 as.integer(n_pos),
                 as.integer(n_center),
                 as.integer(start_d),
                 as.double(marker),
                 coincidence=as.double(rep(0,n_center)),
                 PACKAGE="xoi")

    result <- data.frame(distance=center, coincidence=output$coincidence/n_ind)
    attr(result, "window") <- window
    result[!is.na(result[,2]),]
}
