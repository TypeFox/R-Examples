# Poincare Constant computation for Derivative-based Global Sensitivity Measures (DGSM)           
# Author: Jana Fruth (2014)
#
# Reference: 
# O. Roustant, J. Fruth, B. Iooss and S. Kuhnt,
# Crossed-derivative-based sensitivity measures for interaction screening, 
# Mathematics and Computers in Simulation, 105:105-118, 2014

PoincareConstant <- function(densityfct=dnorm, qfct=qnorm, cdfct, 
                             truncated=FALSE, min=0, max=1, 
                             logconcave=TRUE, optimize.interval=c(-100, 100), ...){
  
  if (logconcave == TRUE){
    if (truncated == FALSE) res <- 1/densityfct(qfct(0.5,...),...)^2
    if (truncated == TRUE){
      res <- (cdfct(max,...)- cdfct(min,...))^2 / (densityfct(qfct((cdfct(min,...)+cdfct(max,...))/2,...),...))^2
    }
  }
  if (logconcave == FALSE){
    fct <- function(x){
      cdf.at.x <- cdfct(x, ...)
      density.at.x <- densityfct(x, ...)
      apply(cbind(cdf.at.x, 1-cdf.at.x),1,min)/(density.at.x)
    }
    c1 <- optimize(f=fct, interval=optimize.interval, maximum=TRUE)$objective
    res <- 4*c1^2
  }
  print(res)
}
