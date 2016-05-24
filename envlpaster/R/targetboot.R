


##########################################################
targetboot <- function(model, nboot, index, u, 
  data, quiet = FALSE, m = 100){
  timer <- proc.time()

  # extract important envelope initial quantities
  # obtain the target and nuisance parameters
  # build U w.r.t. the target parameters
  fam <- model$fam
  pred <- model$pred
  root <- model$root
  x <- model$x
  vars <- colnames(x)
  n <- nrow(x)
  nnode <- ncol(x)
  formula <- model$formula
  modmat.model <- model$modmat
  modelmatrix <- matrix(modmat.model, nrow = n * nnode)
  dimensions <- dim(modmat.model)
  offset <- as.vector(model$origin)

  # obtain tau
  mu <- predict(model, parm.type = "mean.value", 
    model.type = "unconditional")
  tau <- crossprod(modelmatrix, mu)
  p <- length(tau)
  nuis.ind <- c(1:p)[-c(index)]
  target.ind <- index
  k <- length(index)
  nuisance <- tau[nuis.ind]
  target <- tau[target.ind]
  U <- target %o% target

  # run the 1-d algorithm w.r.t. the target
  # construct the envelope estimator
  avar <- model$fisher[target.ind,target.ind]
  foo <- manifold1Dplus(M = avar, U = U, u = u)
  P <- projection(foo)
  tau.env <- crossprod(P, target)
  fulltau <- rep(0,p)
  fulltau[target.ind] <- tau.env
  fulltau[nuis.ind] <- nuisance

  # change the model matrix
  M <- t(modelmatrix)
  M2 <- M[index,]; M2 <- P %*% M2
  M[index,] <- M2
  modelmatrix.int <- t(M)
  modmat.model.int <- array(modelmatrix, dimensions)

  # convert from tau to beta
  beta.foo <- transformUnconditional(parm = fulltau, 
    modelmatrix.int, data, from = "tau", to = "beta", 
    offset = offset, tolerance = 1e-10)

  # set up for the bootstrap for the envelope
  theta.hat <- predict(model, model.type = "cond", 
    parm.type = "canon", newcoef = beta.foo)
  theta.hat <- matrix(theta.hat, nrow = n, ncol = nnode)
  est  <- matrix(nrow = k, ncol = nboot)
  b <- "try-error" 
  class(b) <- "try-error"

  # set up for the bootstrap for the MLE
  theta.hat2 <- predict(model, model.type = "cond", 
    parm.type = "canon")
  theta.hat2 <- matrix(theta.hat2, nrow = n, ncol = nnode)
  est.tau <- matrix(nrow = k, ncol = nboot) # changed from nrow = p
  b2 <- "try-error"
  class(b2) <- "try-error"

  # inital quantities to save computing time
  G.star <- matrix(0, nrow = k, ncol = u)
  P.star <- matrix(0, nrow = k, ncol = k)
  xstar <- xstar2 <- matrix(0, nrow = n, ncol = nnode)
  target.star2 <- tau.env.star <- target.star <- rep(0, k)
  avar.star <- matrix(0, nrow = k, ncol = k)
  modelmatrix.star <- matrix(0, nrow = n*nnode, ncol = p)
  mu.star <- rep(0,nnode*n)
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
    modelmatrix.star <- matrix(aout4star2$modmat, nrow = n*nnode)
    mu.star <- predict(aout4star2, parm.type = "mean.value", 
      model.type = "unconditional")
    est.tau[,iboot] <- crossprod(modelmatrix.star, 
      mu.star)[target.ind]

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
      fam, modmat.model.int, parm = beta.foo), silent = TRUE))[1] 
    
    if(class(b) == "try-error"){
      class(b) <- class(try(aout4star <- aster(xstar, root, pred,
        fam, modmat.model, parm = beta.foo, method = "nlm"), 
        silent = TRUE))[1] 
    }
 
    if(class(b) == "try-error"){
      class(b) <- class(try(aout4star <- aster(xstar, root, pred,
        fam, modmat.model.int), silent = TRUE))[1]
    }

    if(class(b) == "try-error"){
      class(b) <- class(try(aout4star <- aster(xstar, root, pred,
        fam, modmat.model.int, method = "nlm"), silent = TRUE))[1] 
    }

    # get the bootstrapped envelope estimators
    modelmatrix.star <- matrix(aout4star$modmat, nrow = n*nnode)
    mu.star <- predict(aout4star, parm.type = "mean.value", 
      model.type = "unconditional")
    avar.star <- (aout4star$fisher)[target.ind,target.ind]
    G.star <- manifold1Dplus(M = avar.star,
      U = target.star %o% target.star, u = u)
    P.star <- tcrossprod(G.star)
    M <- t(matrix(modmat.model.int, nrow = n*nnode))
    M2 <- M[index,]; M2 <- P.star %*% M2
    M[index,] <- M2
    modelmatrix.star <- t(M)
    tau.env.star <- crossprod(modelmatrix.star, mu.star)
    target.star <- tau.env.star[target.ind]
    est[,iboot] <- crossprod(P.star, target.star)

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
  means.MLE <- apply(est.tau, FUN = mean, 
    MARGIN = 1)
  S <- var(t(est))
  S2 <- var(t(est.tau))
  ratio <- sqrt(diag(S2) / diag(S))
  table <- cbind(tau.env, sqrt(diag(S)), tau[target.ind], 
    sqrt(diag(S2)), ratio)
  colnames(table) <- c("env","se(env)","MLE","se(MLE)","ratio")
  out <- list(u = u, table = table, S = S, S2 = S2,  env.boot.out = est,
    MLE.boot.out = est.tau)
  return(out)
}
##########################################################

















##########################################################
eigenboot <- function(model, nboot, index, vectors,
  data, quiet = FALSE, m = 100){
  timer <- proc.time()

  # extract important envelope initial quantities
  # obtain the target and nuisance parameters
  # build U w.r.t. the target parameters
  u <- length(vectors)
  beta <- model$coef
  fam <- model$fam
  pred <- model$pred
  root <- model$root
  x <- model$x
  vars <- colnames(x)  
  n <- nrow(x)
  nnode <- ncol(x)
  formula <- model$formula
  modmat.model <- model$modmat
  dimensions <- dim(modmat.model)
  modelmatrix <- matrix(modmat.model, nrow = n * nnode)
  offset <- as.vector(model$origin)


  # obtain tau
  mu <- predict(model, parm.type = "mean.value", 
    model.type = "unconditional")
  tau <- crossprod(modelmatrix, mu)
  p <- length(tau)
  nuis.ind <- c(1:p)[-c(index)]
  target.ind <- index
  k <- length(index)
  nuisance <- tau[nuis.ind]
  target <- tau[target.ind]
  U <- target %o% target


  # obtain the eigenspace used to construct the envelope 
  # estimator
  avar <- (model$fisher)[target.ind,target.ind]
  eig <- eigen(avar, symmetric = TRUE)  
  G <- eig$vec[,c(vectors)]
  P <- tcrossprod(G)
  tau.env <- crossprod(P,target)
  fulltau <- rep(0,p)
  fulltau[target.ind] <- tau.env
  fulltau[nuis.ind] <- nuisance


  # change the model matrix
  M <- t(modelmatrix)
  M2 <- M[index,]; M2 <- P %*% M2
  M[index,] <- M2
  modelmatrix.int <- t(M)
  modmat.model.int <- array(modelmatrix, dimensions)


  # convert from tau to beta
  # changed from fulltau to tau
  beta.foo <- transformUnconditional(parm = fulltau, 
    modelmatrix.int, data, from = "tau", to = "beta", 
    offset = offset, tolerance = 1e-20)

  # set up for the bootstrap for the envelope
  theta.hat <- transformUnconditional(parm = fulltau, 
    modelmatrix.int, data, from = "tau", to = "theta", 
    offset = offset, tolerance = 1e-20)
  theta.hat <- matrix(theta.hat, nrow = n, ncol = nnode)
  est  <- matrix(nrow = k, ncol = nboot)
  b <- "try-error"; class(b) <- "try-error"

  # set up for the bootstrap for the MLE
  theta.hat2 <- predict(model, model.type = "cond", 
    parm.type = "canon")
  theta.hat2 <- matrix(theta.hat2, nrow = n, ncol = nnode)
  est.tau <- matrix(nrow = k, ncol = nboot) # changed from nrow = p
  b2 <- "try-error"; class(b2) <- "try-error"

  # inital quantities to save computing time
  G.star <- matrix(0, nrow = k, ncol = u)
  P.star <- matrix(0, nrow = k, ncol = k)
  xstar <- xstar2 <- matrix(0, nrow = n, ncol = nnode)
  tau.env.star <- rep(0, p)
  target.star2 <- target.star <- rep(0, k)
  avar.star <- matrix(0, nrow = k, ncol = k)
  modelmatrix.star <- matrix(0, nrow = n*nnode, ncol = p)
  mu.star <- rep(0,nnode*n)
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
    modelmatrix.star <- matrix(aout4star2$modmat, nrow = n*nnode)
    mu.star <- predict(aout4star2, parm.type = "mean.value", 
      model.type = "unconditional")
    est.tau[,iboot] <- crossprod(modelmatrix.star, 
      mu.star)[target.ind]

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
      fam, modmat.model.int, parm = beta.foo), silent = TRUE))[1] 
    
    if(class(b) == "try-error"){
      class(b) <- class(try(aout4star <- aster(xstar, root, pred,
        fam, modmat.model.int, parm = beta.foo, method = "nlm"), 
        silent = TRUE))[1] 
    }
 
    if(class(b) == "try-error"){
      class(b) <- class(try(aout4star <- aster(xstar, root, pred,
        fam, modmat.model.int), silent = TRUE))[1]
    }

    if(class(b) == "try-error"){
      class(b) <- class(try(aout4star <- aster(xstar, root, pred,
        fam, modmat.model.int, method = "nlm"), silent = TRUE))[1] 
    }

    # get the bootstrapped envelope estimators    
    avar.star <- (aout4star$fisher)[target.ind,target.ind]
    G.star <- eigen(avar.star, symmetric = TRUE)$vec[,c(vectors)]  
    P.star <- tcrossprod(G.star)
    M <- t(matrix(modmat.model.int, nrow = n*nnode))
    M2 <- M[index,]; M2 <- P.star %*% M2
    M[index,] <- M2
    modelmatrix.star <- t(M)

    #modmat.model.star <- array(modelmatrix.star, dimensions)    
    #mu.star <- predict(aout4star, parm.type = "mean.value", 
    #  model.type = "unconditional")
    #tau.env.star <- crossprod(modelmatrix.star, mu.star)
    #target.star <- tau.env.star[target.ind]
    #est[,iboot] <- target.star

    est[,iboot]  <- transformUnconditional(parm = aout4star$coef, 
      modelmatrix.star, data, from = "beta", to = "tau", 
      offset = offset, tolerance = 1e-200)[target.ind]
    

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
  means.MLE <- apply(est.tau, FUN = mean, 
    MARGIN = 1)
  S <- var(t(est))
  S2 <- var(t(est.tau))
  ratio <- sqrt(diag(S2) / diag(S))
  table <- cbind(tau.env, sqrt(diag(S)), tau[target.ind], 
    sqrt(diag(S2)), ratio)
  colnames(table) <- c("env","se(env)","MLE","se(MLE)","ratio")  
  out <- list(u = u, table = table, S = S, S2 = S2, 
    env.boot.out = est, MLE.boot.out = est.tau)
  return(out)
}
##########################################################














