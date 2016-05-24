is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

logit <- function(p) log(p)-log(1-p)
expit <- function(x) exp(x)/(1+exp(x))

lp.lambda.s.y <- function(y, W1, W2, hfunc, obspersubj, pi1.pi0.ratio, gamma.hat, DAT.ods){
  ny <- length(y)
  nsubj <- length(obspersubj)
  if(is.matrix(W1)) nobs <- dim(W1)[1]
  if(is.vector(W1)) nobs <- length(W1)
  z.mat <- matrix(0, nrow=ny*nobs, ncol=(dim(W1)[2] + dim(W2)[2]))
  for(j in 1:ny){
    z.mat[(nobs*(j-1)+1):(nobs*j),] <- cbind(W1, W2*hfunc(y[j]))
  }
  tmpsplit <- rep(0, sum(obspersubj))
  index <- c(0,cumsum(obspersubj))
  for(i in 2:(length(obspersubj)+1)){
    tmpsplit[(index[i-1]+1):index[i]] <- rep(unique(DAT.ods$id)[i-1], each=obspersubj[i-1])
  }
  tmpsplit <- rep(tmpsplit, ny)
  
  firstsplit <- split( t(tcrossprod(gamma.hat,  z.mat)),tmpsplit )
  lpmatrix <- matrix(NA, nrow=nobs, ncol=ny)
  index <- c(1, 1+cumsum(obspersubj[1:(nsubj-1)]))
  for(i in 1:nsubj){
    for(j in 1:obspersubj[i]){
      lpmatrix[index[i]+j-1,] <- split(firstsplit[[i]], rep(1:obspersubj[i], ny))[[j]]
    }
  }
  return(lpmatrix + log(pi1.pi0.ratio))
}



lambda.p.y <- function(y, W1, W2, hfunc, obspersubj, pi1.pi0.ratio, gamma.hat, DAT.ods){
  lp.lambda.s <- lp.lambda.s.y(y, W1, W2, hfunc, obspersubj, pi1.pi0.ratio, gamma.hat, DAT.ods)
  lambda.ok <- lp.lambda.s<710
  ifelse(lambda.ok, expit( apply(lp.lambda.s, 2, "-", DAT.ods$offset.z)) , 1)
}

rho.scaled.y   <- function(lpy, pi1.pi0.ratio){
  1 - lpy + pi1.pi0.ratio*lpy
}


Fy    <- function(lpy, pi1.pi0.ratio){
  Const <- 1- pi1.pi0.ratio
  lpy*(1-lpy)* Const/(1 - lpy + pi1.pi0.ratio*lpy)
}


int <- function(input, points, type){
  stopifnot(is.matrix(input) || dim(input)[2] == length(points))
  stopifnot(is.character(type))
  dorc <- charmatch(type, c("discrete", "continuous"))
  if(dorc==1){
    return(input%*%rep(1, dim(input)[2]))
  }else{
    npt <- length(points)
    a <- points[1:(npt-1)]
    b <- points[2:npt]
    return((input[, 2:npt] + input[,1:(npt-1)])%*%((1/2)*(b-a)))
  }
}


fitz <- function(DAT.ods, w1.formula, w2.formula, y.formula, hfunc){
  nobs <- length(DAT.ods$id)
  terms.w1 <- terms(w1.formula, data=DAT.ods)
  w1.matrix <- model.matrix(w1.formula, data=DAT.ods)[,-1]
  terms.w2 <- terms(w2.formula, data=DAT.ods)
  w2.matrix <- model.matrix(w2.formula, data=DAT.ods)
  terms.y <- terms(y.formula, data=DAT.ods)
  y.matrix <- model.matrix(y.formula, data=DAT.ods)
  
  Y <- model.frame(y.formula, data=DAT.ods)[,attr(terms.y, "response")]
  Z <- model.frame(w1.formula, data=DAT.ods)[,attr(terms.w1, "response")]
  
  for(i in 1:dim(w2.matrix)[2]){attr(w2.matrix, "dimnames")[[2]][i]  <- paste(attr(w2.matrix, "dimnames")[[2]][i], ".w2", sep="")}
  z.index <- attr(terms.w1, "response")
  if(z.index == 0){stop("Must provide a response in the formula for W1")}
  z.df <- data.frame(cbind("z" = Z, w1.matrix, apply(w2.matrix, 2, "*", hfunc(Y))))
  
  z.formula <- paste(names(z.df)[1], "~", names(z.df)[2], sep="")
  count <- 3
  while(count <= length(names(z.df))){
    z.formula <- paste(z.formula, names(z.df)[count], sep="+")
    count <-count+1
  }
  z.df <-data.frame(z.df, "offset.z" = DAT.ods$offset.z, "id" = DAT.ods$id)
  z.formula.offset <- paste(z.formula, "offset(offset.z)", sep="+")
  
  z.formula.offset <- formula(z.formula.offset)
  z.formula        <- formula(z.formula)
  
  optval <- getOption("warn")
  options(warn=-1)
  mod.z            <- glm(z.formula.offset,  data=z.df, family=binomial, x=T)
  options(warn=optval)
  
  W1 <- as.matrix(cbind("intercept" = rep(1, nobs), w1.matrix))
  W2 <- as.matrix(w2.matrix)
  if(!mod.z$converged){stop("Model for referral probability did not converge")}
  return(list(mod.z=mod.z, W1=W1, W2=W2))
}

dmudgamma1 <- function(eta, support, F.y, Fy0, odds.s.y, type, W1, mu.s, hfunc){
  ny <- length(support)
  nobs <- length(Fy0)
  Fy0minusFy <- apply(-F.y, 2, "+", Fy0)
  oddsout <- odds.s.y(support, eta)  
  ymat <- matrix(rep(support, nobs), nrow=nobs, byrow=T)
  
  y.minus.mu <- apply(ymat, 2, "-", mu.s)
  
  y0denom <- int(oddsout, support, type)
  ynum <- as.vector(int(Fy0minusFy*y.minus.mu*oddsout , support, type))
  
  W1timesint <- apply(W1, 2, "*", ynum/y0denom)
  
  as.matrix(W1timesint)
}


dmudgamma2 <- function(eta, support, F.y, Fy0, odds.s.y, type, W2, mu.s, hfunc, y0){
  ny <- length(support)
  nobs <- length(Fy0)
  Fy0minusFy <- apply(-t(apply(F.y, 1, "*", hfunc(support))), 2, "+", hfunc(y0)*Fy0)
  oddsout <- odds.s.y(support, eta)
  
  ymat <- matrix(rep(support, nobs), nrow=nobs, byrow=T)
  
  y.minus.mu <- apply(ymat, 2, "-", mu.s)
  
  y0denom <- int(oddsout, support, type)
  ynum <- as.vector(int(Fy0minusFy*y.minus.mu*oddsout , support, type))
  
  W2timesint <- apply(W2, 2, "*", ynum/y0denom)
  
  as.matrix(W2timesint)
}



getvar <- function(mod.z, mod.y, odds.s.y, VarFun, Y, type, hfunc, support, F.y, Fy0, y0, W1, W2, obspersubj){
  
  mu.s <- fitted(mod.y)
  nobs <- length(mu.s)
  ## Note: scale parameter not needed as it cancels out in sandwich estimate
  sca.param <- mod.y$phi
  
  ### NOW WE MOVE ON TO STANDARD ERRORS ###
  # I calculate as many things with sparse matrices as possible, trying to get all observations
  # accounted for in each matrix
  
  p.z <- length(mod.z$coefficients)
  p.y <- length(mod.y$beta)
  p.all <- p.z+p.y
  
  W  <- mod.z$x
  X  <- mod.y$X
  Z <- mod.z$y
  
  # Y and Z were defined at the beginning
  
  lambda.s <- predict(mod.z, type="response")
  lambda.s[lambda.s==1] <- 1- .Machine$double.eps
  lambda.s[lambda.s==0] <- .Machine$double.eps
  lambda.p    <- ifelse(predict(mod.z) < 709, expit( -1*mod.z$offset.z + predict(mod.z)) , 1)
  
  R.a.inv <- mod.y$R.a.inv*sca.param
  
  bigA <- VarFun(mod.y$eta)
  
  bigT <- apply(W, 2, "*", Z-lambda.s)
  bigU <- t(crossprod( Diagonal(x=sqrt(bigA)) %*%X, R.a.inv %*% Diagonal(x=1/sqrt(bigA)) %*% Diagonal(x=Y-mu.s) ))*sca.param
  
  
  # Calculate derivative of mu w.r.t. gamma2
  
  
  dmdg1 <- dmudgamma1(mod.y$eta, support, F.y, Fy0, odds.s.y, type, W1, mu.s, hfunc)
  dmdg2 <- dmudgamma2(mod.y$eta, support, F.y, Fy0, odds.s.y, type, W2, mu.s, hfunc, y0)
  
  ## Subject specific calculation
  Adiag <- Diagonal(x=sqrt(bigA))
  
  BDiag <- getBlockDiag(obspersubj)$BDiag
  bigTU <- cBind(bigT, bigU)
  Q <- crossprod(bigTU, BDiag) %*% bigTU
  ITTBDiag <- getBlockDiag(rep(p.z, nobs))$BDiag
  
  
  longW <- as.vector(t(W))
  diagones <- matrix(rep(as.vector(Diagonal(x=rep(1, p.z))), nobs), ncol=p.z, byrow=T)
  
  lamdiag <- Diagonal(x = rep(lambda.s*(1-lambda.s), each=p.z))
  ITT <-  crossprod(diagones, Diagonal(x=longW) %*% ITTBDiag %*% lamdiag %*% Diagonal(x=longW) %*% diagones)
  
  IUU <- crossprod(Adiag %*% X, R.a.inv %*% Adiag %*% X)
  
  dgam.mat <- cbind(dmdg1, dmdg2)
  
  IUT <- sca.param*crossprod(Adiag %*% X, R.a.inv %*% solve(Adiag) %*% dgam.mat)
  ITU  <- matrix(0, p.z, p.y)
  I    <- rbind(cbind(as.matrix(ITT), ITU), cbind(as.matrix(IUT), as.matrix(IUU)))  
  
  IQinv <- solve(I,Q)
  
  AVAR <- t(solve(I, t(IQinv)))
  return(sqrt(diag(AVAR)[(p.z+1):(p.z+p.y)])) 
}
