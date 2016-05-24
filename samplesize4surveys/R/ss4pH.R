#' @import TeachingSampling
#' @export
#' 
#' @title
#' The required sample size for testing a null hyphotesis for a single proportion
#' @description 
#' This function returns the minimum sample size required for testing a null hyphotesis regarding a single proportion.
#' @details
#' We assume that it is of interest to test the following set of hyphotesis:
#' \deqn{H_0: P - P_0 = 0 \ \ \ \ vs. \ \ \ \ H_a: P - P_0 = D \neq 0 }
#' Note that the minimun sample size, restricted to the predefined power \eqn{\beta} and confidence \eqn{1-\alpha}, is defined by: 
#' \deqn{n = \frac{S^2}{\frac{D^2}{(z_{1-\alpha} + z_{\beta})^2}+\frac{S^2}{N}}}
#' Where \deqn{S^2=p(1-p)DEFF} 
#'  
#' @author Hugo Andres Gutierrez Rojas <hugogutierrez at usantotomas.edu.co>
#' @param N The population size.
#' @param p The value of the estimated proportion.
#' @param DEFF The design effect of the sample design. By default \code{DEFF = 1}, which corresponds to a simple random sampling design.
#' @param conf The statistical confidence. By default \code{conf = 0.95}.
#' @param power The statistical power. By default \code{power = 0.80}.
#' @param p0 The value to test for the single proportion.
#' @param plot Optionally plot the effect against the sample size.
#' 
#' @references 
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas y estimacion de parametros}. Editorial Universidad Santo Tomas
#' @seealso \code{\link{e4p}}
#' @examples 
#' ss4pH(N = 10000, p = 0.5, p0 = 0.55)
#' ss4pH(N = 10000, p = 0.5, p0 = 0.55, plot=TRUE)
#' ss4pH(N = 10000, p = 0.5, p0 = 0.55, DEFF = 2, plot=TRUE)
#' ss4pH(N = 10000, p = 0.5, p0 = 0.55, conf = 0.99, power = 0.9, DEFF = 2, plot=TRUE)
#' 
#' #############################
#' # Example with BigLucy data #
#' #############################
#' data(BigLucy)
#' attach(BigLucy)
#' 
#' N <- nrow(BigLucy)
#' p <- prop.table(table(SPAM))[1]
#' 
#' # The minimum sample size for testing 
#' # H_0: P - P_0 = 0   vs.   H_a: P - P_0 = D = 0.1
#' D = 0.1 
#' p0 = p - D 
#' ss4pH(N, p, p0, conf = 0.99, power = 0.9, DEFF = 2, plot=TRUE)
#' 
#' # The minimum sample size for testing 
#' # H_0: P - P_0 = 0   vs.   H_a: P - P_0 = D = 0.02
#' D = 0.02
#' p0 = p - D 
#' ss4pH(N, p, p0, conf = 0.99, power = 0.9, DEFF = 3.45, plot=TRUE)


ss4pH = function(N, p, p0, DEFF=1, conf=0.95, power=0.8, plot=FALSE){
  
  S2=p*(1-p)*DEFF
  Za = conf
  Zb = power 
  Z = qnorm(Za)+qnorm(Zb)
  D = abs(p - p0)
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