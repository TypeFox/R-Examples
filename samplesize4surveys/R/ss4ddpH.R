#' @import TeachingSampling
#' @export
#' 
#' @title
#' The required sample size for testing a null hyphotesis for a double difference of proportions
#' @description 
#' This function returns the minimum sample size required for testing a null hyphotesis regarding a double difference of proportion.
#' @details
#' We assume that it is of interest to test the following set of hyphotesis:
#' \deqn{H_0: (P_1 - P_2) - (P_3 - P_4) = 0 \ \ \ \ vs. \ \ \ \ H_a:  (P_1 - P_2) - (P_3 - P_4) = D \neq 0 }
#' Note that the minimun sample size, restricted to the predefined power \eqn{\beta} and confidence \eqn{1-\alpha}, is defined by: 
#' \deqn{n = \frac{S^2}{\frac{D^2}{(z_{1-\alpha} + z_{\beta})^2}+\frac{S^2}{N}}}
#' Where \eqn{S^2=(P_1Q_1+P_2Q_2+P_3Q_3+P_4Q_4)DEFF} and \eqn{Q_i=1-P_i} for \eqn{i=1,2, 3, 4}.
#' @author Hugo Andres Gutierrez Rojas <hugogutierrez at usantotomas.edu.co>
#' @param N The maximun population size between the groups (strata) that we want to compare.
#' @param P1 The value of the first estimated proportion.
#' @param P2 The value of the second estimated proportion.
#' @param P3 The value of the thrid estimated proportion.
#' @param P4 The value of the fourth estimated proportion.
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
#' ss4ddpH(N = 100000, P1 = 0.5, P2 = 0.5, P3 = 0.5, P4 = 0.5, D=0.03)
#' ss4ddpH(N = 100000, P1 = 0.5, P2 = 0.5, P3 = 0.5, P4 = 0.5, D=0.03, plot=TRUE)
#' ss4ddpH(N = 100000, P1 = 0.5, P2 = 0.5, P3 = 0.5, P4 = 0.5, D=0.03, DEFF = 2, plot=TRUE)
#' ss4ddpH(N = 100000, P1 = 0.5, P2 = 0.5, P3 = 0.5, P4 = 0.5, 
#' D=0.03, conf = 0.99, power = 0.9, DEFF = 2, plot=TRUE)
#' 
#' #################################
#' # Example with BigLucyT0T1 data #
#' #################################
#' data(BigLucyT0T1)
#' attach(BigLucyT0T1)
#' 
#' BigLucyT0 <- BigLucyT0T1[Time == 0,]
#' BigLucyT1 <- BigLucyT0T1[Time == 1,]
#' N1 <- table(BigLucyT0$SPAM)[1]
#' N2 <- table(BigLucyT1$SPAM)[1]
#' N <- max(N1,N2)
#' P1 <- prop.table(table(BigLucyT0$ISO))[1]
#' P2 <- prop.table(table(BigLucyT1$ISO))[1]
#' P3 <- prop.table(table(BigLucyT0$ISO))[2]
#' P4 <- prop.table(table(BigLucyT1$ISO))[2]
#' # The minimum sample size for simple random sampling
#' ss4ddpH(N, P1, P2, P3, P4, D = 0.05, plot=TRUE)
#' # The minimum sample size for a complex sampling design
#' ss4ddpH(N, P1, P2, P3, P4, D = 0.05, DEFF = 2, T = 0.5, R = 0.5, conf=0.95, plot=TRUE)

ss4ddpH = function(N, P1, P2, P3, P4, D, DEFF=1, conf=0.95, power=0.8, T = 0, R = 1, plot=FALSE){
  
  Q1 = 1-P1
  Q2 = 1-P2
  Q3 = 1-P3
  Q4 = 1-P4
  S2 <- (P1 * Q1 + P2 * Q2 + P3 * Q3 + P4 * Q4) * (1 - (T * R)) * DEFF
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
      Dseq[k]=100*sqrt(varseq)
    }
    
    plot(nseq,Dseq, type="l", lty=2, pch=1, col=3,ylab="Null effect (D) %",xlab="Sample size")
    points(n.hyp, 100*D, pch=8,bg = "blue")
    abline(h=100*D,lty=3)
    abline(v=n.hyp,lty=3)
  }
  
  result = n.hyp
  result 
}
