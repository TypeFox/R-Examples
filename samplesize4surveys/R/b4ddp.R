#' @import TeachingSampling
#' @export
#' 
#' @title
#' Statistical power for a hyphotesis testing on a difference of proportions
#' @description 
#' This function computes the power for a (right tail) test of difference of proportions. 
#' @return 
#' The power of the test.
#' @details
#' We note that the power is defined as: \deqn{1-\Phi(Z_{1-\alpha} - \frac{(D - [(P_1 - P_2) - (P_3 - P_4)])}{\sqrt{\frac{DEFF}{n}(1-\frac{n}{N})(P_1Q_1+P_2Q_2+P_3Q_3+P_4Q_4)}})}
#' @author Hugo Andres Gutierrez Rojas <hugogutierrez at usantotomas.edu.co>
#' @param N The population size.
#' @param n The sample size.
#' @param P1 The value of the first estimated proportion.
#' @param P2 The value of the second estimated proportion.
#' @param P3 The value of the third estimated proportion.
#' @param P4 The value of the fourth estimated proportion.
#' @param D The value of the null effect.
#' @param DEFF The design effect of the sample design. By default \code{DEFF = 1}, which corresponds to a simple random sampling design.
#' @param conf The statistical confidence. By default \code{conf = 0.95}.
#' @param plot Optionally plot the errors (cve and margin of error) against the sample size.
#' 
#' @references 
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas y estimacion de parametros}. Editorial Universidad Santo Tomas
#' @seealso \code{\link{ss4p}}
#' @examples 
#' b4ddp(N = 10000, n = 400, P1 = 0.5, P2 = 0.5, P3 = 0.5, P4 = 0.5, D = 0.03)
#' b4ddp(N = 10000, n = 400, P1 = 0.5, P2 = 0.5, P3 = 0.5, P4 = 0.5, D = 0.03, plot = TRUE)
#' b4ddp(N = 10000, n = 4000, P1 = 0.5, P2 = 0.5, P3 = 0.5, P4 = 0.5, 
#' D = 0.05, DEFF = 2, conf = 0.99, plot = TRUE)


b4ddp <- function(N, n, P1, P2, P3, P4, D, DEFF = 1, conf = 0.95, plot = FALSE){
  
  Q1 = 1-P1
  Q2 = 1-P2
  Q3 = 1-P3
  Q4 = 1-P4
  S2 = (P1*Q1 + P2*Q2 + P3*Q3 + P4*Q4) * DEFF
  Za = qnorm(conf)
  f = n/N
  VAR = (1 / n) * (1 - f) * S2
  beta = 100 * (1 - pnorm(Za - ((D - ((P1 - P2)-(P3 - P4))) / sqrt(VAR))))
  
  if(plot == TRUE) {
    
    nseq=seq(1,N,10)
    betaseq=rep(NA,length(nseq))
    
    for(k in 1:length(nseq)){
      fseq=nseq[k]/N
      varseq= (1/nseq[k])*(1-fseq)*S2
      betaseq[k]=100*(1 - pnorm(Za - ((D - ((P1 - P2) - (P3 - P4))) / sqrt(varseq))))
    }
    
    plot(nseq,betaseq, type="l", lty=1, pch=1, col=3, ylab="Power of the test (%)",xlab="Sample Size")
    points(n,beta, pch=8, bg = "blue")
    abline(h=beta,lty=3)
    abline(v=n,lty=3)
    
  }
  
  msg <- cat('With the parameters of this function: N =', N, 'n = ', n, 'P1 =', P1, 'P2 =', P2,
             'P3 =', P3, 'P4 =', P4, 'D =', D, 'DEFF = ', DEFF,  'conf =', conf, '. 
             \nThe estimated power of the test is ', beta, '. \n \n')
  
  result <- list(Power = beta)
  result
}
