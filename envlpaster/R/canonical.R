targetbootcanon <- function(model, nboot, index, u, 
  quiet = FALSE, m = 100){
  timer <- proc.time()

  # extract important envelope initial quantities
  # obtain the target and nuisance parameters
  # build U w.r.t. the target parameters
  beta <- model$coef
  p <- length(beta)
  nuis.ind <- c(1:p)[-c(index)]
  target.ind <- index
  k <- length(index)
  nuisance <- beta[nuis.ind]
  target <- beta[target.ind]
  U <- target %o% target

  # important aster model quantities
  formula <- model$formula
  modmat.model <- model$modmat
  fam <- model$fam
  pred <- model$pred
  root <- model$root
  x <- model$x
  n <- nrow(x)
  nnode <- ncol(x)

  # run the 1-d algorithm w.r.t. the target
  # construct the envelope estimator
  avar <- solve(model$fisher, symmetric = TRUE)[target.ind,target.ind]
  foo <- manifold1D.plus(M = avar, U = U, u = u)
  beta.env <- projection(foo) %*% target  
  fullbeta <- rep(0,p)
  fullbeta[target.ind] <- beta.env
  fullbeta[nuis.ind] <- nuisance

  # set up for the bootstrap for the envelope
  theta.hat <- predict(model, model.type = "cond", parm.type = "canon",
    newcoef = fullbeta)
  theta.hat <- matrix(theta.hat, nrow = n, ncol = nnode)
  est  <- matrix(nrow = k, ncol = nboot)
  b <- "try-error" 
  class(b) <- "try-error"

  # set up for the bootstrap for the MLE
  theta.hat2 <- predict(model, model.type = "cond", 
    parm.type = "canon")
  theta.hat2 <- matrix(theta.hat2, nrow = n, ncol = nnode)
  est.beta <- matrix(nrow = k, ncol = nboot) # changed from nrow = p
  b2 <- "try-error"
  class(b2) <- "try-error"

  # inital quantities to save computing time
  G.star <- matrix(0, nrow = k, ncol = u)
  xstar <- xstar2 <- matrix(0, nrow = n, ncol = nnode)
  target.star2 <- beta.env.star <- target.star <- rep(0, k)
  avar.star <- matrix(0, nrow = k, ncol = k)
  aout4star <- aout4star2 <- model


  # the MLE bootstrap
  for(iboot in 1:nboot){

    xstar2 <- raster(theta.hat2, pred, fam, root)
    class(b2) <- class(try(aout4star2 <- aster(xstar2, root, pred,
      fam, modmat.model, parm = beta), silent = TRUE))[1] 
    
    if(class(b2) == "try-error"){
      class(b2) <- class(try(aout4star2 <- aster(xstar2, root, pred,
        fam, modmat.model, parm = beta, 
        method = "nlm"), silent = TRUE))[1] 
    }

    if(class(b2) == "try-error"){
      class(b2) <- class(try( aout4star2 <- aster(xstar2, root, 
        pred, fam, modmat.model), silent = TRUE))[1]
    }

    if(class(b2) == "try-error"){
      class(b2) <- class(try(aout4star2 <- aster(xstar2, root, pred,
        fam, modmat.model, method = "nlm"), silent = TRUE))[1] 
    }

    # get the bootstrapped MLE estimators
    est.beta[,iboot] <- aout4star2$coef[target.ind]

    if(quiet == FALSE){
      if((iboot %% m) == 0){
        timer <- proc.time() - timer
        cat("iteration: ", iboot, " time: ", timer, "\n")
        timer <- proc.time()
      }
    }
  }


  # the envelope bootstrap
  for(iboot in 1:nboot){

    xstar <- raster(theta.hat, pred, fam, root)
    colnames(xstar) <- vars

    # fit the aster model to the regenerated data with the 
    # partial envelope structure imposed
    class(b) <- class(try(aout4star <- aster(xstar, root, pred,
      fam, modmat.model, parm = fullbeta), silent = TRUE))[1] 
    
    if(class(b) == "try-error"){
      class(b) <- class(try(aout4star <- aster(xstar, root, pred,
        fam, modmat.model, parm = fullbeta, method = "nlm"), 
        silent = TRUE))[1] 
    }
 
    if(class(b) == "try-error"){
      class(b) <- class(try(aout4star <- aster(xstar, root, pred,
        fam, modmat.model), silent = TRUE))[1]
    }

    if(class(b) == "try-error"){
      class(b) <- class(try(aout4star <- aster(xstar, root, pred,
        fam, modmat.model, method = "nlm"), silent = TRUE))[1] 
    }

    # get the bootstrapped envelope estimators
    target.star <- aout4star$coef[target.ind]
    avar.star <- solve(aout4star$fisher, symmetric = TRUE,
      tol = 1e-20)[target.ind,target.ind]
    G.star <- manifold1D.plus(M = avar.star, 
      U = target.star %o% target.star, u = u)
    est[,iboot] <- projection(G.star) %*% target.star

    if(quiet == FALSE){
      if((iboot %% m) == 0){
        timer <- proc.time() - timer
        cat("iteration: ", iboot, " time: ", timer, "\n")
        timer <- proc.time()
      }
    }
  }


  # construct the sample variance for the envelope procedure
  # build the output
  means <- apply(est, FUN = mean, MARGIN = 1)
  means.MLE <- apply(est.beta, FUN = mean, 
    MARGIN = 1)
  #means.MLE <- apply(est.beta, FUN = mean, 
  #  MARGIN = 1)[target.ind]
  S <- var(t(est))
  S2 <- var(t(est.beta))
  # S2 <- var.boot[target.ind,target.ind]    
  ratio <- sqrt(diag(S2) / diag(S))
  table <- cbind(beta.env, target, means, means.MLE, ratio)
  colnames(table)[1] <- c("beta.env")
  out <- list(u = u, table = table, S = S, S2 = S2)
  return(out)
}
##########################################################













eigenbootcanon <- function(model, nboot, index, vectors, 
  u, quiet = FALSE, m = 100){
  timer <- proc.time()

  # extract important envelope initial quantities
  # obtain the target and nuisance parameters
  # build U w.r.t. the target parameters
  beta <- model$coef
  p <- length(beta)
  nuis.ind <- c(1:p)[-c(index)]
  target.ind <- index
  k <- length(index)
  nuisance <- beta[nuis.ind]
  target <- beta[target.ind]
  U <- target %o% target

  # important aster model quantities
  formula <- model$formula
  modmat.model <- model$modmat
  fam <- model$fam
  pred <- model$pred
  root <- model$root
  x <- model$x
  n <- nrow(x)
  nnode <- ncol(x)

  # run the 1-d algorithm w.r.t. the target
  # construct the envelope estimator
  avar <- solve(model$fisher, symmetric = TRUE)[target.ind,target.ind]
  eig <- eigen(avar, symmetric = TRUE)  
  G <- eig$vec[,c(vectors)]
  P <- tcrossprod(G)
  beta.env <- crossprod(P,beta[target.ind])
  fullbeta <- rep(0,p)
  fullbeta[target.ind] <- beta.env
  fullbeta[nuis.ind] <- nuisance

  # set up for the bootstrap for the envelope
  theta.hat <- predict(model, model.type = "cond", parm.type = "canon",
    newcoef = fullbeta)
  theta.hat <- matrix(theta.hat, nrow = n, ncol = nnode)
  est  <- matrix(nrow = k, ncol = nboot)
  b <- "try-error" 
  class(b) <- "try-error"

  # set up for the bootstrap for the MLE
  theta.hat2 <- predict(model, model.type = "cond", 
    parm.type = "canon")
  theta.hat2 <- matrix(theta.hat2, nrow = n, ncol = nnode)
  est.beta <- matrix(nrow = k, ncol = nboot) # changed from nrow = p
  b2 <- "try-error"
  class(b2) <- "try-error"

  # inital quantities to save computing time
  xstar <- xstar2 <- matrix(0, nrow = n, ncol = nnode)
  target.star2 <- beta.env.star <- target.star <- rep(0, k)
  avar.star <- matrix(0, nrow = k, ncol = k)
  G.star <- matrix(0, nrow = p, ncol = u)
  P.star <- matrix(0, nrow = p, ncol = p)


  # the MLE bootstrap
  for(iboot in 1:nboot){

    xstar2 <- raster(theta.hat2, pred, fam, root)
    class(b2) <- class(try(aout4star2 <- aster(xstar2, root, pred,
      fam, modmat.model, parm = beta), silent = TRUE))[1] 
    
    if(class(b2) == "try-error"){
      class(b2) <- class(try(aout4star2 <- aster(xstar2, root, pred,
        fam, modmat.model, parm = beta, 
        method = "nlm"), silent = TRUE))[1] 
    }

    if(class(b2) == "try-error"){
      class(b2) <- class(try( aout4star2 <- aster(xstar2, root, 
        pred, fam, modmat.model), silent = TRUE))[1]
    }

    if(class(b2) == "try-error"){
      class(b2) <- class(try(aout4star2 <- aster(xstar2, root, pred,
        fam, modmat.model, method = "nlm"), silent = TRUE))[1] 
    }

    # get the bootstrapped MLE estimators
    est.beta[,iboot] <- aout4star2$coef[target.ind]

    if(quiet == FALSE){
      if((iboot %% m) == 0){
        timer <- proc.time() - timer
        cat("iteration: ", iboot, " time: ", timer, "\n")
        timer <- proc.time()
      }
    }
  }


  # the envelope bootstrap
  for(iboot in 1:nboot){

    xstar <- raster(theta.hat, pred, fam, root)
    colnames(xstar) <- vars

    # fit the aster model to the regenerated data with the 
    # partial envelope structure imposed
    class(b) <- class(try(aout4star <- aster(xstar, root, pred,
      fam, modmat.model, parm = fullbeta), silent = TRUE))[1] 
    
    if(class(b) == "try-error"){
      class(b) <- class(try(aout4star <- aster(xstar, root, pred,
        fam, modmat.model, parm = fullbeta, method = "nlm"), 
        silent = TRUE))[1] 
    }
 
    if(class(b) == "try-error"){
      class(b) <- class(try(aout4star <- aster(xstar, root, pred,
        fam, modmat.model), silent = TRUE))[1]
    }

    if(class(b) == "try-error"){
      class(b) <- class(try(aout4star <- aster(xstar, root, pred,
        fam, modmat.model, method = "nlm"), silent = TRUE))[1] 
    }

    # get the bootstrapped envelope estimators    
    fullbeta.env.star <- aout4star$coef    
    avar.star <- solve(aout4star$fisher, symmetric = TRUE)[target.ind,target.ind]
    eig <- eigen(avar.star, symmetric = TRUE)  
    G.star <- eig$vec[,c(vectors)]
    P.star <- tcrossprod(G.star)
    est[,iboot] <- crossprod(P.star,fullbeta.env.star[target.ind])

    if(quiet == FALSE){
      if((iboot %% m) == 0){
        timer <- proc.time() - timer
        cat("iteration: ", iboot, " time: ", timer, "\n")
        timer <- proc.time()
      }
    }
  }


  # construct the sample variance for the envelope procedure
  # build the output
  means <- apply(est, FUN = mean, MARGIN = 1)
  means.MLE <- apply(est.beta, FUN = mean, 
    MARGIN = 1)
  #means.MLE <- apply(est.beta, FUN = mean, 
  #  MARGIN = 1)[target.ind]
  S <- var(t(est))
  S2 <- var(t(est.beta))
  # S2 <- var.boot[target.ind,target.ind]    
  ratio <- sqrt(diag(S2) / diag(S))
  table <- cbind(beta.env, target, means, means.MLE, ratio)
  colnames(table)[1] <- c("beta.env")
  out <- list(u = u, table = table, S = S, S2 = S2)
  return(out)
}
##########################################################








