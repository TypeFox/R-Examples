#' Convert to SRSF
#'
#' This function converts functions to srsf
#'
#' @param f matrix of functions
#' @param time time
#' @return q matrix of SRSFs
#' @keywords srsf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("simu_data")
#' q = f_to_srvf(simu_data$f,simu_data$time)
f_to_srvf <- function(f,time){
    binsize = mean(diff(time))
    eps = .Machine$double.eps
    tmp = gradient.spline(f,binsize)
    q = tmp$g/sqrt(abs(tmp$g)+eps)
    return(q)
}
