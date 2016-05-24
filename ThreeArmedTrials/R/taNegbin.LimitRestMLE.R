# title Limit of restricted maximum-likelihood estimator in case of negative binomial distributed endpoints
# param rateExp1 A numeric value specifying the rate of the experimental treatment group in the alternative hypothesis
# param rateRef1 A numeric value specifying the rate of the reference treatment group in the alternative hypothesis
# param ratePla1 A numeric value specifying the rate of the placebo treatment group in the alternative hypothesis
# param shape1 A numeric value specifying the shape parameter
# param Delta A numeric value specifying the non-inferiority/superiority margin
# param allocation A (non-empty) vector specifying the sample size allocation (wExp, wRef, wPla)
# return A list containing the following components:
# item{rateExp0, rateRef0, rateRef0}{The limit of the maximum-likelihood estimator for the rates when estimed restricted to the boundary of the null hypothesis}
# item{shape0}{The limit of the maximum-likelihood estimator for the shape parameter when estimed restricted to the boundary of the null hypothesis}
# item{sigma2.rest}{The limit of the maximum-likelihood variance estimator for the Wald-type test when restricted to the boundary of the null hypothesis}
taNegbin.LimitRestMLE <- function(rateExp1, rateRef1, ratePla1, shape1, Delta, allocation = c(1/3, 1/3, 1/3)){
  
  KL.Divergenz <- function(zeta, para.alternat, Delta, w){
    
    rateExp1 <- para.alternat[1]
    rateRef1 <- para.alternat[2]
    ratePla1 <- para.alternat[3]
    shape1 <- para.alternat[4]
    
    xLimit <- 1
    while(pnbinom(q=xLimit, mu=max(ratePla1, rateExp1, rateRef1), size = 1/max(shape1,zeta[3])) < 1){
      xLimit <- xLimit * 5
    }
    
    x <- seq(from=0, to=xLimit)
    
    rateExp <- zeta[1]
    rateRef <- zeta[2]
    shape <- zeta[3]
    if( rateExp <= 0 || rateRef <= 0 || shape <= 0){
      return(Inf)
    }
    ratePla <- ( rateExp - Delta * rateRef ) / (1-Delta) 
    if(ratePla <= 0){
      return(Inf)
    }
    KL.Exp <- sum( ( dnbinom(x, mu = rateExp1, size = 1 / shape1, log = TRUE) - 
                       dnbinom(x, mu = rateExp, size = 1 / shape, log = TRUE) ) * 
                     dnbinom(x, mu = rateExp1, size = 1 / shape1, log = FALSE) )
    KL.Ref <- sum( ( dnbinom(x, mu = rateRef1, size = 1 / shape1, log = TRUE) - 
                       dnbinom(x, mu = rateRef, size = 1 / shape, log = TRUE) ) * 
                     dnbinom(x, mu = rateRef1, size = 1 / shape1, log = FALSE) )
    KL.Pla <- sum( ( dnbinom(x, mu = ratePla1, size = 1 / shape1, log = TRUE) - 
                       dnbinom(x, mu = ratePla, size = 1 / shape, log = TRUE) ) * 
                     dnbinom(x, mu = ratePla1, size = 1 / shape1, log = FALSE) )
    
    return( w[1]*KL.Exp + w[2]*KL.Ref + w[3]*KL.Pla )
  }
  
  opt.KL <- optim(c(rateExp1, rateRef1, shape1), KL.Divergenz, para.alternat = c(rateExp1, rateRef1, ratePla1, shape1), Delta = Delta, w = allocation)$par
  
  rateExp0 <- opt.KL[1]
  rateRef0 <- opt.KL[2]
  ratePla0 <- ( rateExp0 - Delta * rateRef0 ) / ( 1 - Delta )
  shape0 <- opt.KL[3]
  sigma2.rest <- rateExp0 * ( 1 + rateExp0 * shape0 ) / allocation[1] + Delta^2 * rateRef0 * ( 1 + rateRef0 * shape0 ) / allocation[2] + (1-Delta)^2 * ratePla0 * ( 1 + ratePla0 * shape0 ) / allocation[3]
  return(list(
    rateExp0 = rateExp0,
    rateRef0 = rateRef0,
    ratePla0 = ratePla0,
    shape0 = shape0,
    sigma2.rest = sigma2.rest
  )
  )
}
