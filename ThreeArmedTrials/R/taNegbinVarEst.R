# title Variance estimator for Wald-type test in case of negative binomially distributed endpoints
# param xExp A (non-empty) numeric vector of data values coming from the experimental treatment group
# param xRef A (non-empty) numeric vector of data values coming from the reference treatment group
# param xPla A (non-empty) numeric vector of data values coming from the placebo group
# param Delta A numeric value specifying the non-inferiority/superiority margin
# param method A character string determing how the variance for the Wald-type test statistic is estimated, must be \emph{SampleVariance}, or \emph{ML} 
# return A list containing the following components:
# item{rateExp0, rateRef0, rateRef0}{The limit of the maximum-likelihood estimator for the rates when estimed restricted to the boundary of the null hypothesis}
# item{shape0}{The limit of the maximum-likelihood estimator for the shape parameter when estimed restricted to the boundary of the null hypothesis}
# item{sigma2.rest}{The limit of the maximum-likelihood variance estimator for the Wald-type test when restricted to the boundary of the null hypothesis}
taNegbinVarEst <- function(xExp, xRef, xPla, Delta, method = c('RML', 'ML', 'SampleVariance')){
  
  method <- match.arg(method)
  
  # Rate estimators
  rateExp <- mean(xExp)
  rateRef <- mean(xRef)
  ratePla <- mean(xPla)
  
  # Sample size allocation
  nExp <- length(xExp)
  nRef <- length(xRef)
  nPla <- length(xPla)
  n <- nExp + nRef + nPla
  wExp <- nExp / n
  wRef <- nRef / n
  wPla <- nPla / n
  
  if( method == 'SampleVariance'){
    varMmeExp <- var(xExp)
    varMmeRef <- var(xRef)
    varMmePla <- var(xPla)
  }
  if( method == 'ML'){
    shapeMLE <- 1/ .C("newton_Shape", zufallszahlen=as.integer(c(xExp, xRef, xPla)), mean1 = as.double(rateExp), mean2 = as.double(rateRef), mean3 = as.double(ratePla), n1 = as.integer(nExp), n2 = as.integer(nRef), n3 = as.integer(nPla), theta_out = as.double(numeric(1)) )$theta_out
    varMmeExp <- rateExp * (1 + rateExp * shapeMLE)
    varMmeRef <- rateRef * (1 + rateRef * shapeMLE)
    varMmePla <- ratePla * (1 + ratePla * shapeMLE)
  }
  if( method == 'RML'){
    shapeMLE <- 1/ .C("newton_Shape", zufallszahlen=as.integer(c(xExp, xRef, xPla)), mean1 = as.double(rateExp), mean2 = as.double(rateRef), mean3 = as.double(ratePla), n1 = as.integer(nExp), n2 = as.integer(nRef), n3 = as.integer(nPla), theta_out = as.double(numeric(1)) )$theta_out
    
    Lh <- function(theta, xExp, xRef, xPla){
      rateExp <- theta[1]
      rateRef <- theta[2]
      ratePla <- theta[3]
      shape <- theta[4]
      Lh.Exp <- sum( dnbinom(xExp, mu = rateExp, size = 1 / shape, log = TRUE) )
      Lh.Ref <- sum( dnbinom(xRef, mu = rateRef, size = 1 / shape, log = TRUE) )
      Lh.Pla <- sum( dnbinom(xPla, mu = ratePla, size = 1 / shape, log = TRUE) )
      return( -(Lh.Exp + Lh.Ref + Lh.Pla) )
    }
    
    ui <- matrix(c(1,-Delta,(Delta-1),0,diag(1, nrow=4, ncol = 4)), byrow = T, ncol = 4)
    theta <- c(rateExp, Delta*rateExp, rateExp*(1-Delta^2)*Delta/(1-Delta), shapeMLE+0.001)

    COptimPar <- constrOptim( theta=theta , f=Lh, grad = NULL,ui = ui, ci=c(0,0,0,0,0), xExp=xExp, xRef=xRef, xPla=xPla)$par
    rateExpRml <- COptimPar[1]
    rateRefRml <- COptimPar[2]
    ratePlaRml <- COptimPar[3]
    shapeRMLE <- COptimPar[4]

    varMmeExp <- rateExpRml * (1 + rateExpRml * shapeRMLE)
    varMmeRef <- rateRefRml * (1 + rateRefRml * shapeRMLE)
    varMmePla <- ratePlaRml * (1 + ratePlaRml * shapeRMLE)
  }
  
  
  varWaldTest <- varMmeExp / wExp + Delta^2 * varMmeRef / wRef + (1-Delta)^2 * varMmePla / wPla
  return(varWaldTest)
}
