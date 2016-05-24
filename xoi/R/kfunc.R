## kfunc.R

#' estimate Ripley's K function
#'
#' estimate the 1-d version of Ripley's K function
#'
#' @param x list with sorted locations of the data
#' @param d values at which to calculate the function
#' @param lengths lengths of segments studied
#' @param exclude distance to exclude
#' @param tol tolerance value
#'
#' @return data frame with \code{d}, \code{k}, and \code{se}
#'
#' @seealso \code{\link{gammacoi}}, \code{\link{stahlcoi}}, \code{\link{coincidence}}
#' @keywords models
#' @examples
#' L <- 103
#' n <- 2000
#' map1 <- sim.map(L, n.mar=104, anchor=TRUE, include.x=FALSE, eq=TRUE)
#' x <- sim.cross(map1, n.ind=n, m=6, type="bc")
#'
#' xoloc <- find.breaks(x)
#'
#' d <- seq(0, 100, by=0.1)[-1]
#' kf <- kfunc(xoloc, d=d, lengths=rep(L, n))
#'
#' plot(k ~ d, data=kf, type="n", yaxs="i", xaxs="i", las=1,
#'      ylim=c(0, max(kf$k + kf$se)))
#' polygon(c(kf$d, rev(kf$d)), c(kf$k + kf$se, rev(kf$k-kf$se)),
#'         border=NA, col="#add8e650")
#' lines(k ~ d, data=kf)
#'
#' @useDynLib xoi
#' @export
kfunc <-
    function(x, d=seq(0,100,by=0.1), lengths, exclude=0, tol=1e-6)
{
    npt <- sapply(x, length)
    if(!any(npt>0)) stop("Need to have some points.")

    x <- x[npt > 0]
    if(missing(lengths))
        lengths <- sapply(x, max)
    else lengths <- lengths[npt > 0]

    if(length(lengths) != length(x))
        stop("length(lengths) != length(x)")
    if(any(d < exclude)) {
        warning("some d < exclude; these excluded")
        d <- d[d > exclude]
    }

    output <- .C("R_kfunc",
                 as.integer(length(x)),
                 as.integer(sapply(x, length)),
                 as.double(unlist(x)),
                 as.double(lengths),
                 as.integer(length(d)),
                 as.double(d),
                 as.double(exclude),
                 k=as.double(rep(0,length(d))),
                 area=as.double(rep(0,length(d))),
                 rate = as.double(0),
                 as.double(tol),
                 PACKAGE="xoi")

    rate <- output$rate
    k <- output$k
    area = output$area
    result <- data.frame(d=d, k=k, se=1/sqrt(output$area * rate))
    attr(result,"rate") <- rate
    class(result) <- c("kfunc","data.frame")
    result
}
