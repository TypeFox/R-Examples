# function to calculate various expressions for non-linearity in an association computation.
ma.nl <-
function (Y, X)
{
  if (! (is.data.frame(Y) || is.vector(Y)) ) stop(" Y must be a 1 col data.frame or a vector ")
  if (! (is.data.frame(X) || is.vector(X) || is.matrix(X) ) ) stop(" X must be a data.frame ")

  Y <- as.data.frame(Y)
  X <- as.data.frame(X)
  if (ncol(Y) != 1) stop(" Y can only have one variable (column) ")
  if(nrow(Y) != nrow(X)) stop(" Y and X must have the same number of samples (rows)")
  
  cc <- complete.cases(cbind(Y,X))
  Y <- Y[cc,]
  X <- X[cc,]
  
  model <- lm(as.vector(Y)~as.matrix(X))
  # linear association
  Rsq <- summary(model)$r.squared
  # total association linear and non-linear
  A <- ma(cbind(X,Y))$A
  # residual association
  rA <- ma(cbind(X,model$residuals))$A 
  # absolute difference
  AmRsq <- max(0,A-Rsq)
  nl1 <- AmRsq
  # as a preportion of total association
  if(A > 0) nl2 <- AmRsq / A else nl2 <- 0 
  # as a preportion of possible non-linear association
  if(1-Rsq > 0) nl3 <- AmRsq / (1-Rsq) else nl3 <- 0
  return(list(Rsq=Rsq,A=A,rA=rA,nl1=nl1,nl2=nl2,nl3=nl3)) 
}
