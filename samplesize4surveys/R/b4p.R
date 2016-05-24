#' @import TeachingSampling
#' @export
#' 
#' @title
#' Statistical power for a hyphotesis testing on a single proportion
#' @description 
#' This function computes the power for a (right tail) test of proportions. 
#' @return 
#' The power of the test.
#' @details
#' We note that the power is defined as: \deqn{1-\Phi(Z_{1-\alpha} - \frac{(D-P)}{\sqrt{\frac{DEFF}{n}(1-\frac{n}{N})(P (1-P))}})}
#' @author Hugo Andres Gutierrez Rojas <hugogutierrez at usantotomas.edu.co>
#' @param N The population size.
#' @param n The sample size.
#' @param P The value of the first estimated proportion.
#' @param D The value of the null effect. Note that \code{D} must be strictly greater than \code{P}.
#' @param DEFF The design effect of the sample design. By default \code{DEFF = 1}, which corresponds to a simple random sampling design.
#' @param conf The statistical confidence. By default \code{conf = 0.95}.
#' @param plot Optionally plot the errors (cve and margin of error) against the sample size.
#' 
#' @references 
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas y estimacion de parametros}. Editorial Universidad Santo Tomas
#' @seealso \code{\link{ss4p}}
#' @examples 
#' b4p(N = 100000, n = 400, P = 0.5, D = 0.55)
#' b4p(N = 100000, n = 400, P = 0.5, D = 0.9, plot = TRUE)
#' b4p(N = 100000, n = 4000, P = 0.5, D = 0.55, DEFF = 2, conf = 0.99, plot = TRUE)


b4p <- function(N, n, P, D, DEFF = 1, conf = 0.95, plot = FALSE){
  
  Q = 1-P
  S2 = P * Q * DEFF 
  Za = qnorm(conf)
  f = n/N
  VAR = (1 / n) * (1 - f) * S2
  beta = 100 * (1 - pnorm(Za - ((D - P) / sqrt(VAR))))
  
  if(plot == TRUE) {
    
    nseq=seq(1,N,10)
    betaseq=rep(NA,length(nseq))
    
    for(k in 1:length(nseq)){
      fseq=nseq[k]/N
      varseq=(1/nseq[k])*(1-fseq)*S2
      betaseq[k]=100*(1 - pnorm(Za - ( (D - P) / sqrt(varseq))))
    }
    
    plot(nseq,betaseq, type="l", lty=1, pch=1, col=3, ylab="Power of the test (%)",xlab="Sample Size")
    points(n,beta, pch=8, bg = "blue")
    abline(h=beta,lty=3)
    abline(v=n,lty=3)
    
  }
  
  msg <- cat('With the parameters of this function: N =', N, 'n = ', n, 'P =', P,
             'D =', D, 'DEFF = ', DEFF,  'conf =', conf, '. \nThe estimated power of the test is ', beta, '. \n \n')
  
  result <- list(Power = beta)
  result
}
