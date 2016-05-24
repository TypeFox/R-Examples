NULL 
#'
#' Analytic exact solution for One-Dimensional Boussinesq Equation in a two-bounded domain with two constant-value Dirichlet Condition 
#' 
#' @param t time coordinate. 
#' @param x spatial coordinate. Default is \code{seq(from=0,to=L,by=by)}.
#' @param big maximum level of Fourier series considered. Default is 10^7. 
#' @param by see \code{\link{seq}}
#' @param L length of the domain. 
#' @param h1 water surface level at \code{x=0}. Left Dirichlet Bounday Condition.
#' @param h2 water surface level at \code{x=L}. Right Dirichlet Bondary Condition.  
#' @param ks Hydraulic conductivity 
#' @param s drainable pororosity (assumed to be constant) 
#' @param p empirical coefficient to estimate hydraulic diffusivity \eqn{D=ks/(s *(p*h1+(1-p)*h2))}. It ranges between 0 and 1.
#' 
#' @return Solutions for the indicated values of  \code{x} and \code{t}.
#' 
#' @export
#' 
#' @seealso \code{\link{beq.lin.dimensionless}}
#' 
#' @examples 
#' L <- 1000
#' x <- seq(from=0,to=L,by=L/100)
#' t <- 4 # 4 days 
#' h_sol0 <- beq.lin(x=x,t=t*24*3600,h1=2,h2=1,ks=0.01,L=L,s=0.4,big=100,p=0.0)
#' h_solp <- beq.lin(x=x,t=t*24*3600,h1=2,h2=1,ks=0.01,L=L,s=0.4,big=100,p=0.5)
#' h_sol1 <- beq.lin(x=x,t=t*24*3600,h1=2,h2=1,ks=0.01,L=L,s=0.4,big=100,p=1.0)
#' 
#' plot(x,h_sol0,type="l",lty=1,main=paste("Water Surface Elevetion after",t,"days",sep=" "),xlab="x[m]",ylab="h[m]")
#' lines(x,h_solp,lty=2)
#' lines(x,h_sol1,lty=3)
#' legend("topright",lty=1:3,legend=c("p=0","p=0.5","p=1"))

#' 
#' 
#' @author Emanuele Cordano
#' 



beq.lin <-  function(t=0,x=seq(from=0,to=L,by=by),h1=1,h2=1,L=100,ks=0.01,s=0.4,big=10^7,by=L/100,p=0.5) {
  
    D=(ks*(h1*p+h2*(1-p))/s) 
    
    out0 <- h2
    out1 <- (h1-h2)*beq.lin.dimensionless(t=(t*D)/L^2,x=x/L,big=100000)
    
    return(out0+out1)
}


