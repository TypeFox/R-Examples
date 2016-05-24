#' Outlier Detection
#'
#' This function calculates outlier's using geodesic distances of the SRVFs from
#' the median
#'
#' @param q matrix (\eqn{N} x \eqn{M}) of \eqn{M} SRVF functions with \eqn{N}
#' samples
#' @param time vector of size \eqn{N} describing the sample points
#' @param mq median calcuated using \code{\link{time_warping}}
#' @param k cutoff threshold (default = 1.5)
#' @return \item{q_outlier}{outlier functions}
#' @keywords srvf outlier detection
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("toy_data")
#' data("toy_warp")
#' q_outlier = outlier.detection(toy_warp$q0,toy_data$time,toy_warp$mqn,k=.1)
outlier.detection <- function(q, time, mq, k = 1.5){
    N = ncol(q)
    ds = rep(0,N)
    for (kk in 1:N)
        ds[kk] = sqrt(sum(simpson(time, (mq-q[,kk])^2)))

    quartile_range = quantile(ds)
    IQR = quartile_range[4] - quartile_range[2]

    thresh = quartile_range[4] + k * IQR

    ind = which(ds > thresh)

    q_outlier = q[,ind]

    return(q_outlier)
}
