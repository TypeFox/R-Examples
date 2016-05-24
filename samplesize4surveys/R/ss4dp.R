#' @import TeachingSampling
#' @export
#' 
#' @title
#' The required sample size for estimating a single difference of proportions
#' @description 
#' This function returns the minimum sample size required for estimating a single proportion subjecto to predefined errors.
#' @details
#' Note that the minimun sample size to achieve a particular margin of error \eqn{\varepsilon} is defined by: 
#' \deqn{n = \frac{n_0}{1+\frac{n_0}{N}}}
#' Where \deqn{n_0=\frac{z^2_{1-\frac{\alpha}{2}}S^2}{\varepsilon}}
#' and
#' \deqn{S^2=p(1-p)DEFF}
#' Also note that the minimun sample size to achieve a particular coefficient of variation \eqn{cve} is defined by:
#' \deqn{n = \frac{S^2}{p^2cve^2+\frac{S^2}{N}}} 
#'   
#' @author Hugo Andres Gutierrez Rojas <hugogutierrez at usantotomas.edu.co>
#' @param N The maximun population size between the groups (strata) that we want to compare.
#' @param P1 The value of the first estimated proportion.
#' @param P2 The value of the second estimated proportion.
#' @param DEFF The design effect of the sample design. By default \code{DEFF = 1}, which corresponds to a simple random sampling design.
#' @param conf The statistical confidence. By default conf = 0.95. By default \code{conf = 0.95}.
#' @param cve The maximun coeficient of variation that can be allowed for the estimation.
#' @param me The maximun margin of error that can be allowed for the estimation.
#' @param plot Optionally plot the errors (cve and margin of error) against the sample size.
#' 
#' @references 
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas y estimacion de parametros}. Editorial Universidad Santo Tomas
#' @seealso \code{\link{e4p}}
#' @examples 
#' ss4dp(N=100000, P1=0.5, P2=0.55, cve=0.05, me=0.03)
#' ss4dp(N=100000, P1=0.5, P2=0.55, cve=0.05, me=0.03, plot=TRUE)
#' ss4dp(N=100000, P1=0.5, P2=0.55, DEFF=3.45, conf=0.99, cve=0.03, me=0.03, plot=TRUE)
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
#' P1 <- prop.table(table(SPAM))[1]
#' P2 <- prop.table(table(SPAM))[2]
#' # The minimum sample size for simple random sampling
#' ss4dp(N, P1, P2, DEFF=1, conf=0.99, cve=0.03, me=0.03, plot=TRUE)
#' # The minimum sample size for a complex sampling design
#' ss4dp(N, P1, P2, DEFF=3.45, conf=0.99, cve=0.03, me=0.03, plot=TRUE)

ss4dp = function(N, P1, P2, DEFF=1, conf=0.95, cve=0.05, me=0.03, plot=FALSE){
  
  Q1 = 1-P1
  Q2 = 1-P2
  S2 = (P1*Q1 + P2*Q2)*DEFF
  Z = 1-((1-conf)/2)
  n.cve <- S2/((P1-P2)^2*cve^2+(S2/N))
  n0.me <- (qnorm(Z)^2/me^2)*S2
  n.me <- n0.me/(1+(n0.me/N))
  
  if(plot == TRUE) {
    
    nseq=seq(100,N,10)
    cveseq=rep(NA,length(nseq))
    meseq=rep(NA,length(nseq))
    
    for(k in 1:length(nseq)){
      fseq=nseq[k]/N
      varseq=(1/nseq[k])*(1-fseq)*S2
      cveseq[k]=100*sqrt(varseq)/abs(P1-P2)
      meseq[k]=100*qnorm(Z)*sqrt(varseq)
    }
    
    par(mfrow=c(1,2))
    plot(nseq,cveseq, type="l", lty=2, pch=1, col=3,ylab="Coefficient of variation %",xlab="Sample size")
    points(n.cve, 100*cve, pch=8,bg = "blue")
    abline(h=100*cve,lty=3)
    abline(v=n.cve,lty=3)
    
    plot(nseq,meseq, type="l", lty=2, pch=1, col=3,ylab="Margin of error %",xlab="Sample size")
    points(n.me,100*me, pch=8,bg = "red")
    abline(h=100*me,lty=3)
    abline(v=n.me,lty=3)
  }
  
  msg <- cat('With the parameters of this function: N =', N, 'P1 =', P1, 'P2 =', P2, 'DEFF = ',
             DEFF, 'conf =', conf, '.\n
             The estimated sample size to obatin a maximun coefficient of variation of', 100*cve, '% is n=', ceiling(n.cve), '.
             The estimated sample size to obatin a maximun margin of error of', 100*me, '% is n=', ceiling(n.me), '. \n \n')
  
  result <- list(n.cve = ceiling(n.cve), n.me = ceiling(n.me))
  result 
}


