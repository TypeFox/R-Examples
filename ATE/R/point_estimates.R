#We note that we combine the newton raphson for both the simple estimate 
#and the treatment effect of the treated since there is only a small modification
#needed.

#Also note that this function takes care of the case of multiple treatment 
#groups as this is a generic function with an indicator function for given
#treatment. So we can have 2 treatments or k treatments where k is finite.

newton.r<-function (ini, X, Y, Ti, rho, rho1, rho2, FUNu, ATT = FALSE, max.iter = 100, 
                    tol = 0.001, backtrack = TRUE, 
                    bt.a = 0.3, bt.b = 0.5, verbose = TRUE, ...) {
  
  res <- matrix(ini, nrow = 1)
  N <- length(Y)
  umat <- t(apply(X, 1, FUNu))
  if(ATT){
    umat2<- umat[Ti==0,]
    ubar <- apply(umat2, 2, mean)
  }else{
    ubar <- apply(umat, 2, mean)
  }
  dervs<- c(0)
  
  for (i in 2:max.iter) {
    x.val <- res[i - 1, ]
    objectiveValue<- obj(res[i-1,],umat,ubar,Ti,rho,...)
    
    if(objectiveValue <= -1e30){
      warning("The objective function is unbounded, a different rho() function 
              might be needed.")
      lam.hat <- as.vector(tail(res, 1))
      l.t.u <- apply(umat, 1, crossprod, lam.hat)
      nom <- rep(0, N)
      nom[Ti == 1] <- rho1(l.t.u,...)[Ti == 1]/N
      weights <- nom
      return(list("res" = res, "weights" = weights, "conv" = FALSE))
    }
    del.f <- derv.obj(x.val, Ti, rho1, u = umat, ubar = ubar,...)
    H <- derv2.obj(x.val, Ti, rho2, u = umat,...)
    del.x <- -solve(H, del.f)
    nabla.f <- del.f
    dervs<- c(dervs,sum(nabla.f^2))
    
    step <- 1
    if (backtrack) {
      step <- backtrack(bt.a, bt.b, x.val, del.x, nabla.f, 
                        obj, umat, ubar, Ti, rho,...)
    }
    if (verbose) {
      cat(paste("Iteration Number: ", i - 1, "\n"))
      cat(paste("Norm of Derivative: ", sum(nabla.f^2), 
                "\n"))
      cat(paste("Objective value", objectiveValue , "\n"))
    }
    
    res <- rbind(res, x.val + step * del.x)
    
    if (sum(nabla.f^2) < tol) {
      lam.hat <- as.vector(tail(res, 1))
      l.t.u <- apply(umat, 1, crossprod, lam.hat)
      nom <- rep(0, N)
      nom[Ti == 1] <- rho1(l.t.u,...)[Ti == 1]/N
      weights <- nom
      return(list("res" = res, "weights" = weights, "conv" = TRUE))
    }
  }
  
  warning("The algorithm did not converge")
  lam.hat <- as.vector(tail(res, 1))
  l.t.u <- apply(umat, 1, crossprod, lam.hat)
  nom <- rep(0, N)
  nom[Ti == 1] <- rho1(l.t.u,...)[Ti == 1]/N
  weights <- nom
  return(list("res" = res, "weights" = weights, "conv" = FALSE))
}


#########################################################################
#Main function to obtain the simple point estimate
get.est.simple<- function(ini1,ini2, X, Y, Ti, rho, rho1, rho2,
                          FUNu, max.iter = 100,
                          tol = 1e-3, backtrack= TRUE, bt.a = 0.3, bt.b = 0.5,
                          verbose = TRUE, ...){
  
  N<- length(Y)
  if(length(ini1) != length(FUNu(X[1,]))){
    stop("Incorrect length of initial vector")
  }
  if(length(ini2) != length(FUNu(X[1,]))){
    stop("Incorrect length of initial vector")
  }
  
  if(verbose){
    cat("Fitting Newton Raphson for estimating Weights p: \n")
  }
  
  p.hat<- newton.r(ini1, X, Y, Ti, rho, rho1, rho2, FUNu,
                   ATT = FALSE, max.iter= max.iter, tol = tol, 
                   backtrack = backtrack, bt.a = bt.a,
                   bt.b = bt.b,verbose = verbose, ...)
  
  if(verbose){
    cat("\nFitting Newton Raphson for estimating Weights q: \n")
  }
  
  q.hat<- newton.r(ini2, X, Y, 1-Ti, rho, rho1, rho2, FUNu,
                   ATT = FALSE, max.iter = max.iter, tol = tol, 
                   backtrack = backtrack, bt.a = bt.a,
                   bt.b = bt.b, verbose = verbose, ...)
  
  Y_one <-  sum((p.hat$weights*Y)[Ti==1])
  Y_zero<- sum((q.hat$weights*Y)[Ti==0])
  tau.hat<-  Y_one - Y_zero
  
  w.p<- p.hat$weights
  w.q<- q.hat$weights
  
  lam.p<- tail(p.hat$res,1)
  lam.q<- tail(q.hat$res,1)
  
  conv = TRUE
  if(!p.hat$conv | !q.hat$conv){
    warning("Newton Raphson did not converge for atleast one objective function.")
    conv<- FALSE
    #tau.hat<- NaN
  }
  
  res.l<- list("lam.p" = lam.p, "lam.q" = lam.q, "weights.p" = w.p,
               "weights.q" = w.q, "Y1" = Y_one, "Y0"= Y_zero,"tau" = tau.hat,
               "conv" = conv)
  return(res.l)
}

###############################################################################
#Main function to obtain the point estimate for ATT
get.est.ATT<- function(ini2, X, Y, Ti, rho, rho1, rho2,
                       FUNu, max.iter = 100,
                       tol = 1e-3, backtrack= TRUE, bt.a = 0.3, bt.b = 0.5,
                       verbose = TRUE, ...){
  
  N<- length(Y)
  if(length(ini2) != length(FUNu(X[1,]))){
    stop("Incorrect length of initial vector")
  }
  
  if(verbose){
    cat("\nFitting Newton Raphson for estimating Weights q: \n")
  }
  
  q.hat<- newton.r(ini2, X, Y, 1-Ti, rho, rho1, rho2, FUNu,
                   ATT = TRUE, max.iter = max.iter, tol = tol, 
                   backtrack = backtrack, bt.a = bt.a,
                   bt.b = bt.b,verbose = verbose, ...)
  
  Y_one <- mean(Y[Ti==1])
  Y_zero<- sum((q.hat$weights*Y)[Ti==0])
  tau.hat<-  Y_one - Y_zero
  
  w.q<- q.hat$weights
  lam.q<- tail(q.hat$res,1)
  
  conv = TRUE
  if(!q.hat$conv){
    warning("Newton Raphson did not converge for the objective function.")
    conv<- FALSE
    #tau.hat<- NaN
  }
  
  res.l<- list("lam.q" = lam.q, "weights.q" = w.q, 
               "Y1" = Y_one, "Y0"= Y_zero,"tau" = tau.hat,
               "conv" = conv)
  return(res.l)
}

###############################################################################
#Main function to obtain the point estimate for multiple treatment effect
get.est.MT<- function(ini.mat, X, Y, Ti, rho, rho1, rho2,
                       FUNu, max.iter = 100,
                       tol = 1e-3, backtrack= TRUE, bt.a = 0.3, bt.b = 0.5,
                       verbose = TRUE, ...){
  
  N<- length(Y)
  if(ncol(ini.mat) != length(FUNu(X[1,]))){
    stop("Incorrect length of initial vector")
  }
  if(verbose){
    cat("\nFitting Newton Raphson for estimating Weights q: \n")
  }
  
  lam.mat<- matrix(0, ncol = ncol(ini.mat), nrow = length(unique(Ti))  )
  weights.mat<- matrix(0, ncol = N, nrow = length(unique(Ti))  )
  Yj.mat<- numeric(length(unique(Ti)))
  
  for(j in 0:( length(unique(Ti))-1 ) ){
    temp.Ti<- 1*(Ti==j)
    pj.hat<- newton.r(ini.mat[j+1,], X, Y, temp.Ti, rho, rho1, rho2, FUNu,
                      ATT = FALSE, max.iter = max.iter, tol = tol, 
                      backtrack = backtrack, bt.a = bt.a,
                      bt.b = bt.b,verbose = verbose, ...)
    Yj.hat <-  sum((pj.hat$weights*Y)[temp.Ti==1])
    
    Yj.mat[j+1]<- Yj.hat
    lam.mat[j+1,]<- tail(pj.hat$res,1)
    weights.mat[j+1,] <- pj.hat$weights
    
    conv = TRUE
    if(!pj.hat$conv){
      warning("Newton Raphson did not converge for atleast one objective function.")
      conv<- FALSE
      #tau.hat<- NaN
    } 
  }
  
  res.l<- list("lam.mat" = lam.mat, "weights.mat" = weights.mat, 
               "Yj.hat" = Yj.mat, "conv" = conv)
  return(res.l)
}

###############################################################################