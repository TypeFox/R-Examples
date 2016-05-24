#' Calculates two elastic distance
#'
#' This functions calculates the distances between functions,
#' \eqn{D_y} and \eqn{D_x}, where function 1 is aligned to function 2
#'
#' @param f1 sample function 1
#' @param f2 sample function 2
#' @param time sample points of functions
#' @param lambda controls amount of warping (default = 0)
#' @return Returns a list containing \item{Dy}{amplitude distance}
#' \item{Dx}{phase distance}
#' @keywords srvf alignment, distances
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("simu_data")
#' distances = elastic.distance(simu_data$f[,1],simu_data$f[,2],simu_data$time)
elastic.distance <- function(f1,f2,time,lambda = 0){
    q1 = f_to_srvf(f1,time)
    q2 = f_to_srvf(f2,time)
    gam = optimum.reparam(q1,time,q2,time,lambda)
    fw = approx(time,f2,xout=(time[length(time)]-time[1])*gam + time[1])$y
    qw = f_to_srvf(fw,time)
    Dy = sqrt(sum(trapz(time, (qw-q1)^2)))

    M = length(time)
    psi =  sqrt(diff(gam)*(M-1))
    mu = rep(1,M-1)
    Dx  = Re(acos(sum(mu*psi)/(M-1)))
    return(list(Dy=Dy,Dx=Dx))

}
