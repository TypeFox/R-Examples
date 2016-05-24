#' Resample Curve
#'
#' This function resamples a curve to a number of points
#'
#' @param x matrix defining curve (n,T)
#' @param N Number of samples to re-sample curve, N usually is > T
#' @param mode Open ("O") or Closed ("C") curves
#' @return xn matrix defining resampled curve
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("mpeg7")
#' xn = resamplecurve(beta[,,1,1],200)
resamplecurve <- function(x, N=100, mode="O"){
    n = nrow(x)
    T1 = ncol(x)
    xn = matrix(0,n,N)

    delta = rep(0, T1);

    for (r in 2:T1){
        delta[r] = pvecnorm(x[,r]-x[,r-1],2)
    }

    cumdel = cumsum(delta)/sum(delta)
    newdel = seq(0,1,length.out=N)

    for (r in 1:n){
        xn[r,] = spline(cumdel,x[r,],xout=newdel)$y
    }
    
    if (mode=="C"){
      q = curve_to_q(xn)
      qn = project_curve(q)
      xn = q_to_curve(qn)
    }

    return(xn)
}
