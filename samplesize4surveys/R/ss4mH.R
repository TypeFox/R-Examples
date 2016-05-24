#' @import TeachingSampling
#' @export
#' 
#' @title
#' The required sample size for testing a null hyphotesis for a single mean
#' @description 
#' This function returns the minimum sample size required for testing a null hyphotesis regarding a single mean
#' @details
#' We assume that it is of interest to test the following set of hyphotesis:
#' \deqn{H_0: mu - mu_0 = 0 \ \ \ \ vs. \ \ \ \ H_a: mu - mu_0 = D \neq 0 }
#' Note that the minimun sample size, restricted to the predefined power \eqn{\beta} and confidence \eqn{1-\alpha}, is defined by: 
#' \deqn{n = \frac{S^2}{\frac{D^2}{(z_{1-\alpha} + z_{\beta})^2}+\frac{S^2}{N}}}
#' Where \eqn{S^2=\sigma^2 * DEFF} and \eqn{\sigma^2} is the population variance of the varible of interest.
#'  
#' @author Hugo Andres Gutierrez Rojas <hugogutierrez at usantotomas.edu.co>
#' @param N The population size.
#' @param mu The population mean of the variable of interest.
#' @param sigma The population variance of the variable of interest.
#' @param DEFF The design effect of the sample design. By default \code{DEFF = 1}, which corresponds to a simple random sampling design.
#' @param conf The statistical confidence. By default \code{conf = 0.95}.
#' @param power The statistical power. By default \code{power = 0.80}.
#' @param mu0 The value to test for the single mean.
#' @param plot Optionally plot the effect against the sample size.
#' 
#' @references 
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas y estimacion de parametros}. Editorial Universidad Santo Tomas
#' @seealso \code{\link{e4p}}
#' @examples 
#' ss4mH(N = 10000, mu = 500, mu0 = 505, sigma = 100)
#' ss4mH(N = 10000, mu = 500, mu0 = 505, sigma = 100, plot=TRUE)
#' ss4mH(N = 10000, mu = 500, mu0 = 505, sigma = 100, DEFF = 2, plot=TRUE)
#' ss4mH(N = 10000, mu = 500, mu0 = 505, sigma = 100, conf = 0.99, power = 0.9, DEFF = 2, plot=TRUE)
#' 
#' #############################
#' # Example with BigLucy data #
#' #############################
#' data(BigLucy)
#' attach(BigLucy)
#' 
#' N <- nrow(BigLucy)
#' mu <- mean(Income)
#' sigma <- sd(Income)
#' 
#' # The minimum sample size for testing 
#' # H_0: mu - mu_0 = 0   vs.   H_a: mu - mu_0 = D = 15
#' D = 15 
#' mu0 = mu - D 
#' ss4mH(N, mu, mu0, sigma, conf = 0.99, power = 0.9, DEFF = 2, plot=TRUE)
#' 
#' # The minimum sample size for testing 
#' # H_0: mu - mu_0 = 0   vs.   H_a: mu - mu_0 = D = 32
#' D = 32
#' mu0 = mu - D 
#' ss4mH(N, mu, mu0, sigma, conf = 0.99, power = 0.9, DEFF = 3.45, plot=TRUE)


ss4mH = function(N, mu, mu0, sigma, DEFF=1, conf=0.95, power=0.8, plot=FALSE){
  
  S2 <- sigma^2 * DEFF
  Za = conf
  Zb = power 
  Z = qnorm(Za)+qnorm(Zb)
  D = abs(mu - mu0)
  n.hyp = S2/((D^2/Z^2) + (S2/N))
  n.hyp = ceiling(n.hyp)
  
  if(plot == TRUE) {
    
    nseq=seq(100,N,10)
    Dseq=rep(NA,length(nseq))
    
    for(k in 1:length(nseq)){
      fseq=nseq[k]/N
      varseq=(1/nseq[k])*(1-fseq)*S2*(qnorm(Za)+qnorm(Zb))^2
      Dseq[k]=sqrt(varseq)
    }
    
    plot(nseq,Dseq, type="l", lty=2, pch=1, col=3,ylab="Null effect (D)",xlab="Sample size")
    points(n.hyp, D, pch=8,bg = "blue")
    abline(h=D,lty=3)
    abline(v=n.hyp,lty=3)
  }
  
  result = n.hyp
  result 
}