#' @import TeachingSampling
#' @export
#' 
#' @title
#' The required sample size for testing a null hyphotesis for a single difference of proportions
#' @description 
#' This function returns the minimum sample size required for testing a null hyphotesis regarding a single difference of proportions.
#' @details
#' We assume that it is of interest to test the following set of hyphotesis:
#' \deqn{H_0: mu_1 - mu_2 = 0 \ \ \ \ vs. \ \ \ \ H_a: mu_1 - mu_2 = D \neq 0 }
#' Note that the minimun sample size, restricted to the predefined power \eqn{\beta} and confidence \eqn{1-\alpha}, 
#' is defined by: 
#' \deqn{n = \frac{S^2}{\frac{D^2}{(z_{1-\alpha} + z_{\beta})^2}+\frac{S^2}{N}}}
#' where \eqn{S^2=(\sigma_1^2 + \sigma_2^2) * DEFF}
#' @author Hugo Andres Gutierrez Rojas <hugogutierrez at usantotomas.edu.co>
#' @param N The maximun population size between the groups (strata) that we want to compare.
#' @param mu1 The value of the estimated mean of the variable of interes for the first population.
#' @param mu2 The value of the estimated mean of the variable of interes for the second population.
#' @param sigma1 The value of the estimated variance of the variable of interes for the first population.
#' @param sigma2 The value of the estimated mean of a variable of interes for the second population.
#' @param D The minimun effect to test.
#' @param DEFF The design effect of the sample design. By default \code{DEFF = 1}, which corresponds to a simple random sampling design.
#' @param conf The statistical confidence. By default \code{conf = 0.95}.
#' @param power The statistical power. By default \code{power = 0.80}.
#' @param plot Optionally plot the effect against the sample size.
#' 
#' @references 
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas y estimacion de parametros}. Editorial Universidad Santo Tomas
#' @seealso \code{\link{ss4pH}}
#' @examples 
#' ss4dmH(N = 100000, mu1=50, mu2=55, sigma1 = 10, sigma2 = 12, D=3)
#' ss4dmH(N = 100000, mu1=50, mu2=55, sigma1 = 10, sigma2 = 12, D=1, plot=TRUE)
#' ss4dmH(N = 100000, mu1=50, mu2=55, sigma1 = 10, sigma2 = 12, D=0.5, DEFF = 2, plot=TRUE)
#' ss4dmH(N = 100000, mu1=50, mu2=55, sigma1 = 10, sigma2 = 12, D=0.5, DEFF = 2, conf = 0.99, 
#'        power = 0.9, plot=TRUE)
#' 
#' #############################
#' # Example with BigLucy data #
#' #############################
#' data(BigLucy)
#' attach(BigLucy)
#' 
#' N1 <- table(SPAM)[1]
#' N2 <- table(SPAM)[2]
#' N <- max(N1,N2)
#' 
#' BigLucy.yes <- subset(BigLucy, SPAM == "yes")
#' BigLucy.no <- subset(BigLucy, SPAM == "no")
#' mu1 <- mean(BigLucy.yes$Income)
#' mu2 <- mean(BigLucy.no$Income)
#' sigma1 <- sd(BigLucy.yes$Income)
#' sigma2 <- sd(BigLucy.no$Income)
#' 
#' # The minimum sample size for testing 
#' # H_0: mu_1 - mu_2 = 0   vs.   H_a: mu_1 - mu_2 = D = 3
#' D = 3
#' ss4dmH(N, mu1, mu2, sigma1, sigma2, D, DEFF = 2, plot=TRUE)
#' 
#' # The minimum sample size for testing 
#' # H_0: mu_1 - mu_2 = 0   vs.   H_a: mu_1 - mu_2 = D = 3
#' D = 3
#' ss4dmH(N, mu1, mu2, sigma1, sigma2, D, conf = 0.99, power = 0.9, DEFF = 3.45, plot=TRUE)

ss4dmH = function(N, mu1, mu2, sigma1, sigma2, D, DEFF=1, conf=0.95, power=0.8, plot=FALSE){
  
  S2 = (sigma1^2 + sigma2^2) * DEFF
  Za = conf
  Zb = power 
  Z = qnorm(Za)+qnorm(Zb)
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