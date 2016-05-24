NULL 
#'
#' Song et al.'s analytic solution to Boussinesq equation in a 1D semi-infinite domain with a Dirichlet boundary condition
#' 
#' @param t time coordinate. 
#' @param x spatial coordinate. Default is \code{seq(from=0,to=L,by=by)}.
#' @param h1 water surface level or boundary condition coefficient at \code{x=0}. Left Dirichlet Bounday Condition.
#' @param ks Hydraulic conductivity 
#' @param s drainable pororosity (assumed to be constant) 
#' @param nmax order of power series considered for the analytic solution solution. Default is 4. 
#' @param alpha \eqn{\alpha} exponent see Song at al, 2007
#' 
#' 
#' @return The water surface eletion vs time and space obtained by the analytic solution of Boussinesq Equation 
#' 
#' 
#' @author Emanuele Cordano
#' 
#' @note  For major details, see Song at al, 2007
#' 
#' @seealso \code{\link{beq.song.dimensionless}}
#' 
#' @export
#' @references Song, Zhi-yao;Li, Ling;David, Lockington. (2007), "Note on Barenblatt power series solution to Boussinesq equation",Applied Mathematics and Mechanics,
#' \url{http://www.springerlink.com/content/w0u8667772712801/} ,\url{http://dx.doi.org/10.1007/s10483-007-0612-x}
#' 
#' @examples 
#' L <- 1000
#' x <- seq(from=0,to=L,by=L/100)
#' t <- c(4,5,20) #  days 
#' 
#' h_sol1 <- beq.song(t=t[1]*3600*24,x=x,s=0.4,h1=1,ks=0.01,nmax=10,alpha=0)
#' h_sol2 <- beq.song(t=t[2]*3600*24,x=x,s=0.4,h1=1,ks=0.01,nmax=10,alpha=0)
#' h_sol3 <- beq.song(t=t[3]*3600*24,x=x,s=0.4,h1=1,ks=0.01,nmax=10,alpha=0)
#' 	
#' 	
#' plot(x,h_sol1,type="l",lty=1,main="Water Surface Elevetion (Song at's solution) ",xlab="x[m]",ylab="h[m]")
#' lines(x,h_sol2,lty=2)
#' lines(x,h_sol3,lty=3)
#' legend("topright",lty=1:3,legend=paste("t=",t,"days",sep=" "))


#' 
#' 






#
# Author: Emanuele Cordano
#
# Example: Song'and Lockington's analytic solution to boussinesq equation with two-dounded Dirichlet condition
#

beq.song <- function (t=0.5,x=1.0,s=0.4,h1=1,ks=0.01,nmax=4,alpha=1) {
 
    xi <- x*(2*s*(alpha+1)/(h1*ks*t^(alpha+1)))^0.5
    ax <- coefficient.song.solution(n=nmax,lambda=alpha/(alpha+1))
    xi0 <- sum(ax)^(-0.5)
    out <- h1*t^alpha*beq.song.dimensionless(xi=xi,xi0=xi0,a=ax*xi0^2)
 
    return(out)
}
