norm.SOR <- function(y.formula,
                     w1.formula,
                     w2.formula,
                     y0,
                     hfunc, 
                     id,
                     support,
                     pi1.pi0.ratio,  
                     data = parent.frame(),
                     init.beta=NULL,
                     init.sig.2 = 1,
                     est.var = T,
                     CORSTR="independence"){
  #print(pi1.pi0.ratio)
  sigma.2 <- init.sig.2

  DAT.ods <- data
  id <- DAT.ods$id
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
  X <- model.matrix(y.formula,data=DAT.ods)
  
  ## Observed situations where coefficient estimated for Y interactions in model for Z
  ## were very large so that when I calculated linear predictors, they were over 709
  ## which seems to be R's limit for exponentiation
  
  
  lambda.p.y0 <- lambda.p.y(y0, W1, W2, hfunc, obspersubj, pi1.pi0.ratio, gamma.hat, DAT.ods)
  lpy <- lambda.p.y(support, W1, W2, hfunc, obspersubj, pi1.pi0.ratio, gamma.hat, DAT.ods)
  
  rho.scaled.y0   <- as.vector(1-lambda.p.y0+pi1.pi0.ratio*lambda.p.y0)
  rho.y <- rho.scaled.y(lpy, pi1.pi0.ratio)
  
  ## Calculate F_i (y) for y=1 and y=0
  Fy0    <- Fy(lambda.p.y0, pi1.pi0.ratio)
  F.y <- Fy(lpy, pi1.pi0.ratio)
  
  

  
  InvLink <- function(eta){
    oddsout <- odds.s.y(support, eta)
    yoddsout <- t(apply(odds.s.y(support, eta), 1, "*", support))
    intyodds <- int(yoddsout, support, "cont")
    intodds <- int(oddsout, support, "cont")
    as.vector(intyodds/intodds)
  }
  
  LinkFun <- identity
  
  
  VarFun <- function(eta){
    oddsout <- odds.s.y(support, eta)
    y2oddsout <- t(apply(odds.s.y(support, eta), 1, "*", support^2))
    yoddsout <- t(apply(odds.s.y(support, eta), 1, "*", support))
    inty2odds <- int(y2oddsout, support, "cont")
    intyodds <- int(yoddsout, support, "cont")
    intodds <- int(oddsout, support, "cont")
    as.vector(inty2odds/intodds - (intyodds/intodds)^2)
  }
  
  InvLinkDeriv <- function(eta){
    VarFun(eta)
  }
  
  FunList <- list(LinkFun, VarFun, InvLink, InvLinkDeriv)
  
  y.minus.mu.2 <- function(sigma.2, mu.s){
    rhomat <- apply(log( rho.y ), 2, "-", log( rho.scaled.y0))
    ymat.2 <- apply(matrix(rep(support, nobs), nrow=nobs, byrow=T), 2, "-", mu.s)^2
    toint <- ymat.2*(exp(-ymat.2/(2*sigma.2) + rhomat))
    int(toint, support, "cont")
  }
  
  y.minus.mu.4 <- function(sigma.2, mu.s){
    rhomat <- apply(log( rho.y ), 2, "-", log( rho.scaled.y0))
    ymat.2 <- apply(matrix(rep(support, nobs), nrow=nobs, byrow=T), 2, "-", mu.s)^2
    ymat.4 <- apply(matrix(rep(support, nobs), nrow=nobs, byrow=T), 2, "-", mu.s)^4
    toint <- (ymat.4)*(exp(-ymat.2/(2*sigma.2)+rhomat))
    int(toint, support, "cont")
  }
  
  normalize <- function(sigma.2, mu.s){
    rhomat <- apply(log( rho.y ), 2, "-", log( rho.scaled.y0))
    ymat.2 <- apply(matrix(rep(support, nobs), nrow=nobs, byrow=T), 2, "-", mu.s)^2
    toint <- (exp(-ymat.2/(2*sigma.2)+rhomat))
    int(toint, support, "cont")
  }
  
  if(est.var){
    f.sig <- function(sigma.2, mu.p){
      sum((Y-mu.p)^2) - sum(y.minus.mu.2(sigma.2, mu.p)/normalize(sigma.2, mu.p))
    }
    
    df.sig <- function(sigma.2, mu.p){
      -(1/(2*sigma.2^2))*(sum(y.minus.mu.4(sigma.2, mu.p)/normalize(sigma.2, mu.p)) + (1/(2*sigma.2^2))*sum((y.minus.mu.2(sigma.2, mu.p)^2)/(normalize(sigma.2, mu.p)^2)))
      
    }
    sigstop <- F
    sigma.2.new <- 0
    
    
    
    
    while(!sigstop){
      odds.s.y <- function(y, eta){
        rhomat <- apply(log( rho.y ), 2, "-", log( rho.scaled.y0))
        ymat <- apply(matrix(rep(y, nobs), nrow=nobs, byrow=T), 2, "-", y0)
        etaymat <- apply(ymat, 2, "*", eta)
        referencedist <- -(1/2)*(matrix(rep(y^2, nobs), nrow=nobs, byrow=T)) + (1/2)*(y0^2)
        exp( etaymat/sigma.2 + referencedist/sigma.2 + rhomat )
      }
      
      
      mod.y <- geemR(y.formula, family=FunList, id = id, data=DAT.ods, corstr = CORSTR, Mv = 1, init.beta=init.beta, tol=1e-6, scale.fix=T, init.phi=sigma.2)
      
      #### SOLVE NEWTON RAPHSON
      
      init.beta <- mod.y$beta
      mu.p <- X%*%init.beta
      
      sigma.2.new <- sigma.2 - f.sig(sigma.2, mu.p)/df.sig(sigma.2, mu.p)
      
      if(abs(sigma.2.new - sigma.2)/sigma.2 < 0.001){
        sigstop <- T
      }
      sigma.2 <- sigma.2.new
    }
  }
  odds.s.y <- function(y, eta){
    rhomat <- apply(log( rho.y ), 2, "-", log( rho.scaled.y0))
    ymat <- apply(matrix(rep(y, nobs), nrow=nobs, byrow=T), 2, "-", y0)
    etaymat <- apply(ymat, 2, "*", eta)
    referencedist <- -(1/2)*(matrix(rep(y^2, nobs), nrow=nobs, byrow=T)) + (1/2)*(y0^2)
    exp( etaymat/sigma.2 + referencedist/sigma.2 + rhomat )
  }
  
  
  
  
  mod.y <- geemR(y.formula, family=FunList, id = id, data=DAT.ods, corstr = CORSTR, Mv = 1, init.beta=init.beta, tol=1e-6, scale.fix=T, init.phi=sigma.2)
  
  se.coefs <- getvar(mod.z, mod.y, odds.s.y, VarFun, Y,"c", hfunc, support, F.y, Fy0, y0, W1, W2, obspersubj)
  
  return(list("coefs" = mod.y$beta, "se.coefs" = se.coefs))
}
