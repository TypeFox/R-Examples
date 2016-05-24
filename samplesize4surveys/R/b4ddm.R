#' @import TeachingSampling
#' @export
#' 
#' @title
#' Statistical power for a hyphotesis testing on a double difference of means.
#' @description 
#' This function computes the power for a (right tail) test of double difference of means 
#' @return 
#' The power of the test.
#' @details
#' We note that the power is defined as: \deqn{1-\Phi(Z_{1-\alpha} - \frac{(D - [(\mu_1 - \mu_2) - (\mu_3 - \mu_4)])}{\sqrt{\frac{1}{n}(1-\frac{n}{N})S^2}})}
#' where \deqn{S^2 = DEFF (\sigma_1^2 + \sigma_2^2 + \sigma_3^2 + \sigma_4^2}
#' @author Hugo Andres Gutierrez Rojas <hugogutierrez at usantotomas.edu.co>
#' @param N The population size.
#' @param n The sample size.
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
#' @param D The value of the null effect.
#' @param DEFF The design effect of the sample design. By default \code{DEFF = 1}, which corresponds to a simple random sampling design.
#' @param conf The statistical confidence. By default \code{conf = 0.95}.
#' @param plot Optionally plot the errors (cve and margin of error) against the sample size.
#' 
#' @references 
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas y estimacion de parametros}. Editorial Universidad Santo Tomas
#' @seealso \code{\link{ss4p}}
#' @examples 
#' b4ddm(N = 100000, n = 400, mu1=50, mu2=55, mu3=50, mu4=55, 
#' sigma1 = 10, sigma2 = 12, sigma3 = 10, sigma4 = 12, D = 7)
#' b4ddm(N = 100000, n = 400, mu1=50, mu2=55, mu3=50, mu4=65, 
#' sigma1 = 10, sigma2 = 12, sigma3 = 10, sigma4 = 12, D = 12, plot = TRUE)
#' b4ddm(N = 100000, n = 4000, mu1=50, mu2=55, mu3=50, mu4=65, 
#' sigma1 = 10, sigma2 = 12, sigma3 = 10, sigma4 = 12, D = 11, DEFF = 2, conf = 0.99, plot = TRUE)


b4ddm <- function(N, n, mu1, mu2, mu3, mu4, sigma1, sigma2, sigma3, sigma4, D, DEFF = 1, conf = 0.95, T = 0, R = 1, plot = FALSE){
  S2 <- (sigma1^2 + sigma2^2 + sigma3^2 + sigma4^2) * (1 - (T * R)) * DEFF
  Za = qnorm(conf)
  f = n/N
  VAR = (1 / n) * (1 - f) * S2
  beta = 100 * (1 - pnorm(Za - ((D - ((mu1 - mu2)-(mu3 - mu4))) / sqrt(VAR))))
  
  if(plot == TRUE) {
    
    nseq=seq(1,N,10)
    betaseq=rep(NA,length(nseq))
    
    for(k in 1:length(nseq)){
      fseq=nseq[k]/N
      varseq=(1/nseq[k])*(1-fseq)*S2
      betaseq[k]=100*(1 - pnorm(Za - ((D - ((mu1 - mu2)-(mu3 - mu4))) / sqrt(varseq))))
    }
    
    plot(nseq,betaseq, type="l", lty=1, pch=1, col=3, ylab="Power of the test (%)",xlab="Sample Size")
    points(n,beta, pch=8, bg = "blue")
    abline(h=beta,lty=3)
    abline(v=n,lty=3)
    
  }
  
  msg <- cat('With the parameters of this function: N =', N, 'n = ', n, 'mu1 =',"mu1 =", mu1, "mu2 =", mu2,
             "mu3 =", mu3, "mu4 =", mu4, "sigma1 =", sigma1, "sigma2 =", sigma2, "sigma3 =", sigma3, "sigma4 =", 
             'D =', D, 'DEFF = ', DEFF,  'conf =', conf, 
             '. \nThe estimated power of the test is ', beta, '. \n \n')
  
  result <- list(Power = beta)
  result  
}
