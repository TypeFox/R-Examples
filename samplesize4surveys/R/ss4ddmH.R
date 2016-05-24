#' @import TeachingSampling
#' @export
#' 
#' @title
#' The required sample size for testing a null hyphotesis for a double difference of proportions
#' @description 
#' This function returns the minimum sample size required for testing a null hyphotesis regarding a double difference of proportions.
#' @details
#' We assume that it is of interest to test the following set of hyphotesis:
#' \deqn{H_0: (mu_1 - mu_2) - (mu_3 - mu_4) = 0 \ \ \ \ vs. \ \ \ \ H_a: (mu_1 - mu_2) - (mu_3 - mu_4) = D \neq 0 }
#' Note that the minimun sample size, restricted to the predefined power \eqn{\beta} and confidence \eqn{1-\alpha}, 
#' is defined by: 
#' \deqn{n = \frac{S^2}{\frac{D^2}{(z_{1-\alpha} + z_{\beta})^2}+\frac{S^2}{N}}}
#' where \eqn{S^2=(\sigma_1^2 + \sigma_2^2 + \sigma_3^2 + \sigma_4^2) * DEFF}
#' @author Hugo Andres Gutierrez Rojas <hugogutierrez at usantotomas.edu.co>
#' @param N The maximun population size between the groups (strata) that we want to compare.
#' @param mu1 The value of the estimated mean of the variable of interes for the first population.
#' @param mu2 The value of the estimated mean of the variable of interes for the second population.
#' @param mu3 The value of the estimated mean of the variable of interes for the third population.
#' @param mu4 The value of the estimated mean of the variable of interes for the fourth population.
#' @param sigma1 The value of the estimated variance of the variable of interes for the first population.
#' @param sigma2 The value of the estimated mean of a variable of interes for the second population.
#' @param sigma3 The value of the estimated variance of the variable of interes for the third population.
#' @param sigma4 The value of the estimated mean of a variable of interes for the fourth population.
#' @param T The overlap between waves. By default \code{T = 0}.
#' @param R The correlation between waves. By default \code{R = 1}.
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
#' ss4ddmH(N = 100000, mu1=50, mu2=55, mu3=50, mu4=65, 
#' sigma1 = 10, sigma2 = 12, sigma3 = 10, sigma4 = 12, D=3)
#' ss4ddmH(N = 100000, mu1=50, mu2=55, mu3=50, mu4=65, 
#' sigma1 = 10, sigma2 = 12, sigma3 = 10, sigma4 = 12, D=1, plot=TRUE)
#' ss4ddmH(N = 100000, mu1=50, mu2=55, mu3=50, mu4=65, 
#' sigma1 = 10, sigma2 = 12, sigma3 = 10, sigma4 = 12, D=0.5, DEFF = 2, plot=TRUE)
#' ss4ddmH(N = 100000, mu1=50, mu2=55, mu3=50, mu4=65, 
#' sigma1 = 10, sigma2 = 12, sigma3 = 10, sigma4 = 12, D=0.5, DEFF = 2, conf = 0.99, 
#'        power = 0.9, plot=TRUE)
#' 
#' #############################
#' # Example with BigLucy data #
#' #############################
#' data(BigLucyT0T1)
#' attach(BigLucyT0T1)
#' 
#' BigLucyT0 <- BigLucyT0T1[Time == 0,]
#' BigLucyT1 <- BigLucyT0T1[Time == 1,]
#' N1 <- table(BigLucyT0$ISO)[1]
#' N2 <- table(BigLucyT0$ISO)[2]
#' N <- max(N1,N2)
#' 
#' BigLucyT0.yes <- subset(BigLucyT0, ISO == "yes")
#' BigLucyT0.no <- subset(BigLucyT0, ISO == "no")
#' BigLucyT1.yes <- subset(BigLucyT1, ISO == "yes")
#' BigLucyT1.no <- subset(BigLucyT1, ISO == "no")
#' mu1 <- mean(BigLucyT0.yes$Income)
#' mu2 <- mean(BigLucyT0.no$Income)
#' mu3 <- mean(BigLucyT1.yes$Income)
#' mu4 <- mean(BigLucyT1.no$Income)
#' sigma1 <- sd(BigLucyT0.yes$Income)
#' sigma2 <- sd(BigLucyT0.no$Income)
#' sigma3 <- sd(BigLucyT1.yes$Income)
#' sigma4 <- sd(BigLucyT1.no$Income)
#' 
#' # The minimum sample size for testing 
#' # H_0: (mu_1 - mu_2) - (mu_3 - mu_4) = 0   vs.   
#' # H_a: (mu_1 - mu_2) - (mu_3 - mu_4) = D = 3
#'
#' ss4ddmH(N, mu1, mu2, mu3, mu4, sigma1, sigma2, sigma3, sigma4,
#'  D = 3, conf = 0.99, power = 0.9, DEFF = 3.45, plot=TRUE)

ss4ddmH = function(N, mu1, mu2, mu3, mu4, sigma1, sigma2, sigma3, sigma4, D, DEFF=1, conf=0.95, power=0.8, T = 0, R = 1, plot=FALSE){
  
  S2 = (sigma1^2 + sigma2^2 + sigma3^2 + sigma4^2) * (1 - (T * R)) * DEFF
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
    
    plot(nseq,Dseq, type="l", lty=2, pch=1, col=3,ylab="Null effect",xlab="Sample size")
    points(n.hyp, D, pch=8,bg = "blue")
    abline(h=D,lty=3)
    abline(v=n.hyp,lty=3)
  }
  
  result = n.hyp
  result 
}