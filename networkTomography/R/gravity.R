
#' Gravity estimation for a single time point
#' 
#' @param yt length-m numeric vector of observed aggregate flows at time t
#' @param srcDstInd list of source and destination flow indices corresponding to
#'   each point-to-point flow, as produced by \code{\link{getSrcDstIndices}}
#' @return xhat, a numeric vector of length k providing gravity estimates of the
#'   point-to-point flows of interest
#' @keywords models multivariate ts
#' @export
#' @family gravity
#' @examples
#' data(cmu)
#' srcDstInd <- getSrcDstIndices(cmu$A.full)
#' estimate <- gravity.fit(yt=cmu$Y.full[1,], srcDstInd=srcDstInd)
gravity.fit <- function(yt, srcDstInd) {
    total <- (sum(yt[unique(srcDstInd$src)]) +
              sum(yt[unique(srcDstInd$dst)])) / 2
    p.src <- yt[srcDstInd$src] / total
    p.dst <- yt[srcDstInd$dst] / total
    xhat <- p.src * p.dst * total
    return(xhat)
}

#' Run tomogravity estimation on complete time series of aggregate flows
#' 
#' @param Y n x m matrix contain one vector of observed aggregate flows per row
#' @param srcDstInd list of source and destination flow indices corresponding to
#'   each point-to-point flow, as produced by \code{\link{getSrcDstIndices}}
#' @return Xhat, a n x k matrix containing a vector of estimated point-to-point
#'   flows (for each time point) per row
#' @keywords models multivariate ts
#' @export
#' @family gravity
#' @examples
#' data(cmu)
#' srcDstInd <- getSrcDstIndices(cmu$A.full)
#' estimate <- gravity(Y=cmu$Y[1:3,], srcDstInd=srcDstInd)
gravity <- function(Y, srcDstInd) {
    return(t(apply(Y, 1, gravity.fit, srcDstInd=srcDstInd)))
}
