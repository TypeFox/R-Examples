NULL 
#'
#' Alogoritm for resolution of the series coefficient \eqn{a_n} for the dimensionless formula for \eqn{H} in \code{\link{beq.song.dimensionless}}
#' 
#' @param n approcimation order 
#' @param lambda dimensionless parameter releted to \eqn{\alpha} see Song at al, 2007 
#' 
#' 
#' @return the \eqn{a_n} series coefficient 
#' 
#' @note For major details, see Song at al, 2007
#' 
#' @references Song, Zhi-yao;Li, Ling;David, Lockington. (2007), "Note on Barenblatt power series solution to Boussinesq equation",Applied Mathematics and Mechanics,
#' \url{http://www.springerlink.com/content/w0u8667772712801/} ,\url{http://dx.doi.org/10.1007/s10483-007-0612-x}
#' 
#' @export
#' @author Emanuele Cordano
#' 
#' 
#' 




coefficient.song.solution <- function(n=4,lambda=0) {
  a <- array(NA,n)
  a[1]=1/4
  a[2]=(2*lambda-1)/16
  for (i in 3:n) {
    a[i]=(2*lambda+1-i)/i^2*a[i-1]
    for (k in 2:(i-1)) a[i]=a[i]-2*(i+1)/i*a[k]*a[i+1-k]
 
    
  } 
  return(a)
}  

