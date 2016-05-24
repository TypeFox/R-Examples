
pois.SOR <- function(y.formula,
                      w1.formula,
                      w2.formula,
                      y0,
                      hfunc, 
                      id,
                      support,
                      pi1.pi0.ratio,  
                      data = parent.frame(),
                      init.beta=NULL,
                      CORSTR="independence"){

  # Process the arguments
  DAT.ods <- data 
  ny <- length(support)  
  DAT.ods <- data
  id <- DAT.ods$id
  nobs <- length(id)
  nsubj <- length(unique(id))
  obspersubj <- as.numeric(summary(split(id, id))[,1])
  maxtime <- max(obspersubj)
  pi1.pi0.ratio <- DAT.ods$pi.ratio
  
  aux <- fitz(DAT.ods, w1.formula, w2.formula, y.formula, hfunc)
  
  mod.z <- aux$mod.z
  W1 <- aux$W1
  W2 <- aux$W2
  gamma.hat <- mod.z$coefficients
  
  terms.y <- terms(y.formula, data=DAT.ods)
  y.matrix <- model.matrix(y.formula, data=DAT.ods)
  
  Y <- model.frame(y.formula, data=DAT.ods)[,attr(terms.y, "response")]
  Z <- mod.z$y
  
  
  # Calculate linear predictor for \lambda_S
  # Argument y is a vector, in almost all cases the same as support
  
  
  lambda.p.y0 <- lambda.p.y(y0, W1, W2, hfunc, obspersubj, pi1.pi0.ratio, gamma.hat, DAT.ods)
  lpy <- lambda.p.y(support, W1, W2, hfunc, obspersubj, pi1.pi0.ratio, gamma.hat, DAT.ods)
  
  rho.scaled.y0   <- as.vector(1-lambda.p.y0+pi1.pi0.ratio*lambda.p.y0)
  rho.y <- rho.scaled.y(lpy, pi1.pi0.ratio)
  
  ## Calculate F_i (y) for y=1 and y=0
  Fy0    <- Fy(lambda.p.y0, pi1.pi0.ratio)
  F.y <- Fy(lpy, pi1.pi0.ratio)
  
  
  # The odds funcion, includes the poisson stuff as reference distribution
  odds.s.y <- function(y, eta){
    ny <- length(y)
    rhomat <- log( apply(rho.y ,2,"/", rho.scaled.y0))
    ymat <- apply(matrix(rep(y, nobs), nrow=nobs, byrow=T), 2, "-", y0)
    etaymat <- apply(ymat, 2, "*", eta)
    referencedist <- -lfactorial(matrix(rep(y, nobs), nrow=nobs, byrow=T)) + lfactorial(y0)
    exp( etaymat + referencedist + rhomat )
  }
  
  
  # Inverse link is int(y*odds)/int(odds)
  InvLink <- function(eta){
    oddsout <- odds.s.y(support, eta)
    as.vector((oddsout%*%support)/(oddsout%*%rep(1, ny)))
  }
  
  LinkFun <- log
  
  # Variance function is int((y^2)*odds)/int(odds) - (int(y*odds)/int(odds))^2
  VarFun <- function(eta){
    oddsout <- odds.s.y(support, eta)
    as.vector((oddsout%*%(support^2))/(oddsout%*%rep(1, ny))) - as.vector((oddsout%*%support)/(oddsout%*%rep(1, ny)))^2
  }
  
  # InvLinkDeriv is the same as Variance function
  InvLinkDeriv <- function(eta){
    VarFun(eta)
  }
  
  FunList <- list(LinkFun, VarFun, InvLink, InvLinkDeriv)
  
  # Run through geemR to calculate coefficients
  mod.y <- geemR(y.formula, family=FunList, id = id, data=DAT.ods, corstr = CORSTR, Mv = 1, init.beta = init.beta, scale.fix=T, init.phi=1)
  
  se.coefs <- getvar(mod.z, mod.y, odds.s.y, VarFun, Y,"d", hfunc, support, F.y, Fy0, y0, W1, W2, obspersubj)
  
  return(list("coefs" = mod.y$beta, "se.coefs" = se.coefs))
}