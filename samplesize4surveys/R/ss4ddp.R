#' @import TeachingSampling
#' @export
#' 
#' @title
#' The required sample size for estimating a double difference of proportions
#' @description 
#' This function returns the minimum sample size required for estimating a double difference of proportion subjecto to predefined errors.
#' @details
#' Note that the minimun sample size (for each group at each wave) to achieve a particular margin of error \eqn{\varepsilon} is defined by: 
#' \deqn{n = \frac{n_0}{1+\frac{n_0}{N}}}
#' Where \deqn{n_0=\frac{z^2_{1-\frac{\alpha}{2}}S^2}{\varepsilon}}
#' and
#' \deqn{S^2 = (P1 * Q1 + P2 * Q2 + P3 * Q3 + P4 * Q4) * (1 - (T * R)) * DEFF}
#' Also note that the minimun sample size to achieve a particular coefficient of variation \eqn{cve} is defined by:
#' \deqn{n = \frac{S^2}{(ddp)^2cve^2+\frac{S^2}{N}}} 
#' And \eqn{ddp} is the expected estimate of the double difference of proportions.   
#' @author Hugo Andres Gutierrez Rojas <hugogutierrez at usantotomas.edu.co>
#' @param N The population size.
#' @param P1 The value of the first estimated proportion at first wave.
#' @param P2 The value of the second estimated proportion at first wave.
#' @param P3 The value of the first estimated proportion at second wave.
#' @param P4 The value of the second estimated proportion at second wave.
#' @param DEFF The design effect of the sample design. By default \code{DEFF = 1}, which corresponds to a simple random sampling design.
#' @param T The overlap between waves. By default \code{T = 0}.
#' @param R The correlation between waves. By default \code{R = 1}.
#' @param conf The statistical confidence. By default conf = 0.95. By default \code{conf = 0.95}.
#' @param cve The maximun coeficient of variation that can be allowed for the estimation.
#' @param me The maximun margin of error that can be allowed for the estimation.
#' @param plot Optionally plot the errors (cve and margin of error) against the sample size.
#' 
#' @references 
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas y estimacion de parametros}. Editorial Universidad Santo Tomas
#' @seealso \code{\link{ss4dp}}
#' @examples 
#' ss4ddp(N=100000, P1=0.05, P2=0.55, P3= 0.5, P4= 0.6, cve=0.05, me=0.03)
#' ss4ddp(N=100000, P1=0.05, P2=0.55, P3= 0.5, P4= 0.6, cve=0.05, me=0.03, plot=TRUE)
#' ss4ddp(N=100000, P1=0.05, P2=0.55, P3= 0.5, P4= 0.6, DEFF=3.45, conf=0.99, 
#' cve=0.03, me=0.03, plot=TRUE)
#' ss4ddp(N=100000, P1=0.05, P2=0.55, P3= 0.5, P4= 0.6, DEFF=3.45, conf=0.99,
#'  cve=0.03, me=0.03, T = 0.5, R = 0.9, plot=TRUE)
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
#' ss4ddp(N, P1, P2, P3, P4, conf=0.95, cve=0.05, me=0.03, plot=TRUE)
#' # The minimum sample size for a complex sampling design
#' ss4ddp(N, P1, P2, P3, P4, T = 0.5, R = 0.5, conf=0.95, cve=0.05, me=0.03, plot=TRUE)


ss4ddp<-function(N, P1, P2, P3 ,P4, DEFF = 1, conf = 0.95, cve = 0.05, me = 0.03, T = 0, R = 1, plot = FALSE) 
{
  Q1 <- 1 - P1
  Q2 <- 1 - P2
  Q3 <- 1 - P3
  Q4 <- 1 - P4
  
  S2    <- (P1 * Q1 + P2 * Q2 + P3 * Q3 + P4 * Q4) * (1 - (T * R)) * DEFF
  Z     <- 1 - ((1 - conf)/2)
  n0.me <- (me^2)/(qnorm(Z)^2)
  n.me  <- S2/((n0.me)+(S2/N)) 

  n.cve <- S2/(((P1 - P2) - (P3 - P4))^2 * cve^2 + (S2/N))
  
  if (plot == TRUE & all(diff(c(P1,P2,P3,P4))!=0)) {
    nseq = seq(100, N, 10)
    cveseq = rep(NA, length(nseq))
    meseq = rep(NA, length(nseq))
    for (k in 1:length(nseq)) {
      fseq = nseq[k]/N
      varseq = (1/nseq[k]) * (1 - fseq) * S2
      cveseq[k] = 100 * sqrt(varseq)/abs(P1 - P2)
      meseq[k] = 100 * qnorm(Z) * sqrt(varseq)
    }
    par(mfrow = c(1, 2))
    plot(nseq, cveseq, type = "l", lty = 2, pch = 1, col = 3, 
         ylab = "Coefficient of variation %", xlab = "Sample size")
    points(n.cve, 100 * cve, pch = 8, bg = "blue")
    abline(h = 100 * cve, lty = 3)
    abline(v = n.cve, lty = 3)
    plot(nseq, meseq, type = "l", lty = 2, pch = 1, col = 3, 
         ylab = "Margin of error %", xlab = "Sample size")
    points(n.me, 100 * me, pch = 8, bg = "red")
    abline(h = 100 * me, lty = 3)
    abline(v = n.me, lty = 3)
  }
  msg <- cat("With the parameters of this function: N =", N, 
             "P1 =", P1, "P2 =", P2,"P3 =",P3,"P4 =",P4, "DEFF = ", DEFF, "conf =", conf,"T =",T,"R =",R, 
             ".\n\n The estimated sample size (for each group at each wave) to obatin a maximun coefficient of variation of", 
             100 * cve, "% is n=", ceiling(n.cve), ".\n The estimated sample size (for each group at each wave) to obatin a maximun margin of error of", 
             100 * me, "% is n=", ceiling(n.me), ". \n \n")
  result <- list(n.cve = ceiling(n.cve), n.me = ceiling(n.me))
  result
  
  if(all(diff(c(P1,P2,P3,P4))==0)){
  cat("n.cve and plot are missing because P1, P2, P3 and P4 are all the same")
  }
}

