#' @import TeachingSampling
#' @export
#' 
#' @title
#' The required sample size for estimating a single difference of proportions
#' @description 
#' This function returns the minimum sample size required for estimating a single proportion subjecto to predefined errors.
#' @details
#' Note that the minimun sample size to achieve a relative margin of error \eqn{\varepsilon} is defined by: 
#' \deqn{n = \frac{n_0}{1+\frac{n_0}{N}}}
#' Where \deqn{n_0=\frac{z^2_{1-\frac{alpha}{2}}S^2}{\varepsilon^2 (\mu_1 - \mu_2)^2}} and 
#' \eqn{S^2=(\sigma_1^2 + \sigma_2^2) * DEFF}
#' Also note that the minimun sample size to achieve a coefficient of variation \eqn{cve} is defined by:
#' \deqn{n = \frac{S^2}{|\bar{y}_1-\bar{y}_2|^2 cve^2 + \frac{S^2}{N}}} 
#'   
#' @author Hugo Andres Gutierrez Rojas <hugogutierrez at usantotomas.edu.co>
#' @param N The maximun population size between the groups (strata) that we want to compare.
#' @param mu1 The value of the estimated mean of the variable of interes for the first population.
#' @param mu2 The value of the estimated mean of the variable of interes for the second population.
#' @param sigma1 The value of the estimated variance of the variable of interes for the first population.
#' @param sigma2 The value of the estimated mean of a variable of interes for the second population.
#' @param DEFF The design effect of the sample design. By default \code{DEFF = 1}, which corresponds to a simple random sampling design.
#' @param conf The statistical confidence. By default conf = 0.95. By default \code{conf = 0.95}.
#' @param cve The maximun coeficient of variation that can be allowed for the estimation.
#' @param rme The maximun relative margin of error that can be allowed for the estimation.
#' @param plot Optionally plot the errors (cve and margin of error) against the sample size.
#' 
#' @references 
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas y estimacion de parametros}. Editorial Universidad Santo Tomas
#' @seealso \code{\link{e4p}}
#' @examples 
#' ss4dm(N=100000, mu1=50, mu2=55, sigma1 = 10, sigma2 = 12, cve=0.05, rme=0.03)
#' ss4dm(N=100000, mu1=50, mu2=55, sigma1 = 10, sigma2 = 12, cve=0.05, rme=0.03, plot=TRUE)
#' ss4dm(N=100000, mu1=50, mu2=55, sigma1 = 10, sigma2 = 12, DEFF=3.45, conf=0.99, cve=0.03, 
#'      rme=0.03, plot=TRUE)
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
#' # The minimum sample size for simple random sampling
#' ss4dm(N, mu1, mu2, sigma1, sigma2, DEFF=1, conf=0.99, cve=0.03, rme=0.03, plot=TRUE)
#' # The minimum sample size for a complex sampling design
#' ss4dm(N, mu1, mu2, sigma1, sigma2, DEFF=3.45, conf=0.99, cve=0.03, rme=0.03, plot=TRUE)

ss4dm = function(N, mu1, mu2, sigma1, sigma2, DEFF=1, conf=0.95, cve=0.05, rme=0.03, plot=FALSE){
  
  S2 = (sigma1^2 + sigma2^2)*DEFF
  Z = 1-((1-conf)/2)
  n.cve <- S2 / ((mu1 - mu2)^2 * cve^2 + (S2 / N))
  me <- rme * abs(mu1 - mu2)
  n0 <- (qnorm(Z)^2/me^2)*S2
  n.rme <- n0/(1+(n0/N))
  
  if(plot == TRUE) {
    
    nseq=seq(100,N,10)
    cveseq=rep(NA,length(nseq))
    meseq=rep(NA,length(nseq))
    
    for(k in 1:length(nseq)){
      fseq=nseq[k]/N
      varseq=(1/nseq[k])*(1-fseq)*S2
      cveseq[k]=100*sqrt(varseq)/abs(mu1-mu2)
      meseq[k]=100*qnorm(Z)*sqrt(varseq)
    }
    
    par(mfrow=c(1,2))
    plot(nseq,cveseq, type="l", lty=2, pch=1, col=3,ylab="Coefficient of variation %",xlab="Sample size")
    points(n.cve, 100*cve, pch=8,bg = "blue")
    abline(h=100*cve,lty=3)
    abline(v=n.cve,lty=3)
    
    plot(nseq,meseq, type="l", lty=2, pch=1, col=3,ylab="Margin of error %",xlab="Sample size")
    points(n.rme,100*me, pch=8,bg = "red")
    abline(h=100*me,lty=3)
    abline(v=n.rme,lty=3)
  }
  
  msg <- cat('With the parameters of this function: N =', N, 'mu1 =', mu1, 'mu2 =', mu2, 
             'sigma1 = ', sigma1, 'sigma2 = ', sigma2, 'DEFF = ', DEFF, 'conf =', conf, '.\n
             The estimated sample size to obatin a maximun coefficient of variation of', 100*cve, '% is n=', ceiling(n.cve), '.
             The estimated sample size to obatin a maximun margin of error of', 100*rme, '% is n=', ceiling(n.rme), '. \n \n')
  
  result <- list(n.cve = ceiling(n.cve), n.rme = ceiling(n.rme))
  result 
}



