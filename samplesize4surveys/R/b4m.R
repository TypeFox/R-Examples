#' @import TeachingSampling
#' @export
#' 
#' @title
#' Statistical power for a hyphotesis testing on a single mean
#' @description 
#' This function computes the power for a (right tail) test of means. 
#' @return 
#' The power of the test.
#' @details
#' We note that the power is defined as: \deqn{1-\Phi(Z_{1-\alpha} - \frac{(D - \mu)}{\sqrt{\frac{1}{n}(1-\frac{n}{N})S^2}})}
#' where \deqn{S^2 = DEFF \sigma^2}
#' @author Hugo Andres Gutierrez Rojas <hugogutierrez at usantotomas.edu.co>
#' @param N The population size.
#' @param n The sample size.
#' @param mu The value of the estimated mean of the variable of interest.
#' @param sigma The value of the standard deviation of the variable of interest.
#' @param D The value of the null effect. Note that \code{D} must be strictly greater than \code{mu}.
#' @param DEFF The design effect of the sample design. By default \code{DEFF = 1}, which corresponds to a simple random sampling design.
#' @param conf The statistical confidence. By default \code{conf = 0.95}.
#' @param plot Optionally plot the errors (cve and margin of error) against the sample size.
#' 
#' @references 
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas y estimacion de parametros}. Editorial Universidad Santo Tomas
#' @seealso \code{\link{ss4p}}
#' @examples 
#' b4m(N = 100000, n = 400, mu = 3, sigma = 1, D = 3.1)
#' b4m(N = 100000, n = 400, mu = 5, sigma = 10, D = 7, plot = TRUE)
#' b4m(N = 100000, n = 400, mu = 50, sigma = 100, D = 100, DEFF = 3.4, conf = 0.99, plot = TRUE)


b4m <- function(N, n, mu, sigma, D, DEFF = 1, conf = 0.95, plot = FALSE){
  
  S2 = sigma^2 * DEFF 
  Za = qnorm(conf)
  f = n/N
  VAR = (1 / n) * (1 - f) * S2
  beta = 100 * (1 - pnorm(Za - ( (D - mu) / sqrt(VAR))))
  
  if(plot == TRUE) {
    
    nseq=seq(1,N,10)
    betaseq=rep(NA,length(nseq))
    
    for(k in 1:length(nseq)){
      fseq=nseq[k]/N
      varseq=(1/nseq[k])*(1-fseq)*S2
      betaseq[k]=100*(1 - pnorm(Za - ((D - mu) / sqrt(varseq))))
    }
    
    plot(nseq,betaseq, type="l", lty=1, pch=1, col=3, ylab="Power of the test (%)",xlab="Sample Size")
    points(n,beta, pch=8, bg = "blue")
    abline(h=beta,lty=3)
    abline(v=n,lty=3)
    
  }
  
  msg <- cat('With the parameters of this function: N =', N, 'n = ', n, 'mu =', mu, 'sigma =', sigma,
             'D =', D, 'DEFF = ', DEFF,  'conf =', conf, '. \nThe estimated power of the test is ', beta, '. \n \n')
  
  result <- list(Power = beta)
  result
}
