# kwak_intensity.R

#' Estimate intensity function
#'
#' Estimate intensity function for a chromosome.
#'
#'
#' @param cross Cross object; must be a backcross.  See
#' \code{\link[qtl]{read.cross}} for format details.
#' @param chr Chromosome to consider (only one is allowed).  If missing, the
#' first chromosome is considered.
#' @param window Window size
#' @param ncalc Total number of points for calculations.
#' @return Data frame with columns \code{position} and \code{intensity}.  The
#' input argument \code{window} is kept as an attribute.
#' @author Il youp Kwak
#' @seealso \code{\link{coincidence}}
#' @keywords utilities
#' @examples
#'
#' map1 <- sim.map(103, n.mar=104, anchor=TRUE, include.x=FALSE, eq=TRUE)
#' x <- sim.cross(map1, n.ind=2000, m=6, type="bc")
#'
#' out <- intensity(x)
#' plot(out, type="l", lwd=2, ylim=c(0, max(out[,2])))
#'
#' @useDynLib xoi
#' @export
intensity <-
    function(cross, chr, window=2.5, ncalc=500)
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

    intensitySUB(xomat, window, map, n.ind, ncalc)
}


intensitySUB <-
    function(xomat, window, marker, n_ind, N)
{

    n_xo <- ncol(xomat)
    n_pos <- length(marker)

    center <- seq(marker[1], marker[n_pos], length=N)

    n_center <- length(center)
    xovec <- c(xomat)

    output <- .C("R_get_intensity",
                 as.integer(xovec),
                 as.double(window),
                 as.double(center),
                 as.integer(n_pos),
                 as.integer(n_xo),
                 as.integer(n_center),
                 as.double(marker),
                 intensity=as.double(rep(0,n_center)),
                 as.integer(n_ind),
                 PACKAGE="xoi")

    result <- data.frame(position=center,
                         intensity=output$intensity)
    attr(result, "window") <- window
    result
}
