NULL 
#'
#' Dimensionless solution for one-dimensional derived equation from scaling Boussinesq Equation (Song et al, 2007)
#' 
#' 
#' @param xi dimensionless coordinate (see \code{Note})
#' @param xi0 displacement of wetting front expressed as dimensionless coordinate (see \code{Note})
#' @param a vector of coefficient returned by \code{\link{coefficient.song.solution}}
#' 
#' 
#' @note The expession for the dimensionless coordinate (Song at al., 2007) is \eqn{  \xi=x (\frac{2 \, s }{\eta_1 \, K_s \, t^{\alpha+1} } )^{1/2}}
#' and the solution for the dimensionless equation derived by Boussinesq Equation is: 
#' \eqn{H = \sum_{n=0}^{\infty} a_n (1-\frac{\xi}{\xi_0} )^n} for \eqn{ \xi<\xi_0 } , otherwise is 0 . 
#' 
#' @author Emanuele Cordano
#' 
#' @return the dimesioneless solution, i.e. the variable \eqn{H}
#'
#' @seealso \code{\link{beq.song}}
#' @export 
#' @references Song, Zhi-yao;Li, Ling;David, Lockington. (2007), "Note on Barenblatt power series solution to Boussinesq equation",Applied Mathematics and Mechanics,
#' \url{http://www.springerlink.com/content/w0u8667772712801/} ,\url{http://dx.doi.org/10.1007/s10483-007-0612-x}
 

#
# Author: Emanuele Cordano
#
# Example: Song'and Lockington's analytic solution to boussinesq equation with two-dounded Dirichlet condition
#

beq.song.dimensionless <- function(xi,xi0,a) {
  

  temp <- 1-xi/xi0
  temp[temp<0] <- 0
  
  a0 <- 0
  out <- array(a0,length(temp))
  
  for (n in 1:length(a)) out <- out+a[n]*temp^n
  return(out)
}
