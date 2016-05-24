#' @import TeachingSampling
#' @export
#' 
#' @title
#' The required sample size for estimating a double difference of means
#' @description 
#' This function returns the minimum sample size required for estimating a double difference of means subjecto to predefined errors.
#' @details
#' Note that the minimun sample size to achieve a relative margin of error \eqn{\varepsilon} is defined by: 
#' \deqn{n = \frac{n_0}{1+\frac{n_0}{N}}}
#' Where \deqn{n_0=\frac{z^2_{1-\frac{alpha}{2}}S^2}{\varepsilon^2 \mu^2}} and 
#' \eqn{S^2=(\sigma_1^2 + \sigma_2^2 + \sigma_3^2 + \sigma_4^2) * DEFF}
#' Also note that the minimun sample size to achieve a coefficient of variation \eqn{cve} is defined by:
#' \deqn{n = \frac{S^2}{|(\bar{y}_1-\bar{y}_2) - (\bar{y}_3-\bar{y}_4) |^2 cve^2 + \frac{S^2}{N}}} 
#'   
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
#' ss4ddm(N=100000, mu1=50, mu2=55, mu3=50, mu4=65, 
#' sigma1 = 10, sigma2 = 12, sigma3 = 10, sigma4 = 12, cve=0.05, rme=0.03)
#' ss4ddm(N=100000, mu1=50, mu2=55, mu3=50, mu4=65, 
#' sigma1 = 10, sigma2 = 12, sigma3 = 10, sigma4 = 12, cve=0.05, rme=0.03, plot=TRUE)
#' ss4ddm(N=100000, mu1=50, mu2=55, mu3=50, mu4=65, 
#' sigma1 = 10, sigma2 = 12, sigma3 = 10, sigma4 = 12, DEFF=3.45, conf=0.99, cve=0.03, 
#'      rme=0.03, plot=TRUE)
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
#' # The minimum sample size for simple random sampling
#' ss4ddm(N, mu1, mu2, mu3, mu4, sigma1, sigma2, sigma3, sigma4, 
#' DEFF=1, conf=0.95, cve=0.001, rme=0.001, plot=TRUE)
#' # The minimum sample size for a complex sampling design
#' ss4ddm(N, mu1, mu2, mu3, mu4, sigma1, sigma2, sigma3, sigma4, 
#' DEFF=3.45, conf=0.99, cve=0.03, rme=0.03, plot=TRUE)

ss4ddm = function(N, mu1, mu2, mu3, mu4, sigma1, sigma2, sigma3, sigma4, DEFF=1, conf=0.95, cve=0.05, rme=0.03, T = 0, R = 1, plot=FALSE){
  
  S2 = (sigma1^2 + sigma2^2 + sigma3^2 + sigma4^2) * (1 - (T * R)) * DEFF
  Z = 1-((1-conf)/2)
  n.cve <- S2 / (((mu1 - mu2) - (mu3 - mu4))^2 * cve^2 + (S2 / N))
  me <- rme * abs((mu1 - mu2) - (mu3 - mu4))
  n0 <- (qnorm(Z)^2/me^2)*S2
  n.rme <- n0/(1+(n0/N))
  
  if(plot == TRUE) {
    
    nseq=seq(100,N,10)
    cveseq=rep(NA,length(nseq))
    meseq=rep(NA,length(nseq))
    
    for(k in 1:length(nseq)){
      fseq=nseq[k]/N
      varseq=(1/nseq[k])*(1-fseq)*S2
      cveseq[k]=100*sqrt(varseq)/abs((mu1 - mu2) - (mu3 - mu4))
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
  
  msg <- cat('With the parameters of this function: N =', N, 'mu1 =', mu1, 'mu2 =', mu2, 'mu3 =', mu3, 'mu4 =', mu4,  
             'sigma1 = ', sigma1, 'sigma2 = ', sigma2, 'sigma3 = ', sigma3, 'sigma4 = ', sigma4, 'DEFF = ', DEFF, 'conf =', conf, '.\n
             The estimated sample size to obatin a maximun coefficient of variation of', 100*cve, '% is n=', ceiling(n.cve), '.
             The estimated sample size to obatin a maximun margin of error of', 100*rme, '% is n=', ceiling(n.rme), '. \n \n')
  
  result <- list(n.cve = ceiling(n.cve), n.rme = ceiling(n.rme))
  result 
}



