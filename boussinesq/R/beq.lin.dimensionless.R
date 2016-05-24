NULL 
#'
#' Analytic exact solution for Dimentionless (i. e. diffusivity equal to 1 - unity) One Dimensional Heat Equation in a two-bounded domain with two constant-value Dirichlet Conditions 
#' 
#' @param t time coordinate. 
#' @param x spatial coordinate. Default is \code{seq(from=0,to=L,by=by)}.
#' @param big maximum level of Fourier series considered. Default is 100000. 
#' @param by see \code{\link{seq}}
#' @param L length of the domain.  It is used if \code{x} is not specified.
#' 
#' @return  Solutions for the specificied values of  \code{x} and \code{t}
#' 
#' @references Rozier-Cannon, J. (1984), The One-Dimensional Heat Equation, Addison-Wesley Publishing Company, Manlo Park, California, encyclopedia of Mathematics and its applications.
#' 
#' @seealso \code{\link{beq.lin}}
#' 
#' @export
#' @author Emanuele Cordano

beq.lin.dimensionless <- function(t=0,x=seq(from=0,to=L,by=by),big=100000,by=L*0.01,L=1) {
#
# t is dimensionless and resacaled with L^2/D  and q*L^2/
# x is dimensionless and resacaled with L
# coefficient is dimnsionless and rescaled with (q*L^2)/(D*s)
#

 sum=0
 for (n in 1:big) sum=sum-2/(pi*n)*exp(-n^2*pi^2*t)*sin(n*pi*x) 
 #sum=sum+1/((2*n+1)*pi)^3*(1-exp(-(2*n+1)^2*pi^2*t))*sin((2*n+1)*x*pi)

  sum=sum+1-x
  
  return(sum)

}
