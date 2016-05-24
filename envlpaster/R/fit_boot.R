

fit.boot <- function(model, nboot, index, vectors = NULL, u = NULL,
    data, amat, newdata, modmat.new = NULL, renewdata = NULL, 
    fit.name = NULL, method = c("eigen","1d"), quiet = FALSE, m = 100){

  # stopping condition
  stopifnot(!(is.null(vectors) & is.null(u)))
  stopifnot(!(is.null(fit.name) & is.null(renewdata)))
  timer <- proc.time()

  # extract important envelope initial quantities
  # obtain the target and nuisance parameters
  # build U w.r.t. the target parameters
  beta <- model$coef
  x <- model$x
  vars <- colnames(x)  
  code <- data$code
  families <- data$families
  fam <- model$fam
  pred <- model$pred
  root <- model$root  
  n <- nrow(x)
  nnode <- ncol(x)
  npop <- nrow(newdata)
  p <- length(beta)
  formula <- model$formula
  modmat.model <- model$modmat
  dimensions <- dim(modmat.model)
  modelmatrix <- matrix(modmat.model, nrow = n * nnode)
  offset <- as.vector(model$origin)

  
  # obtain tau
  mu <- predict(model, parm.type = "mean.value", 
    model.type = "unconditional")
  tau <- crossprod(modelmatrix, mu)
  nuis.ind <- c(1:p)[-c(index)]
  target.ind <- index
  k <- length(index)
  nuisance <- tau[nuis.ind]
  target <- tau[target.ind]


  # obtain the projection used to construct the envelope 
  # estimator
  avar <- (model$fisher)[target.ind,target.ind]
  U <- P <- NULL
  if(method == "eigen"){
    u <- length(vectors)
    eig <- eigen(avar, symmetric = TRUE)  
    G <- eig$vec[,c(vectors)]
    P <- tcrossprod(G)  
  }
  if(method == "1d"){
    U <- target %o% target
    G <- manifold1Dplus(M = avar, U = U, u = u)
    P <- tcrossprod(G)
  }

  tau.env <- crossprod(P,target)
  fulltau <- rep(0,p)
  fulltau[target.ind] <- tau.env
  fulltau[nuis.ind] <- nuisance


  # change the model matrix
  M <- t(modelmatrix)
  M2 <- P %*% M[index,]
  M[index,] <- M2
  modelmatrix.int <- t(M)
  modmat.model.int <- array(modelmatrix, dimensions)


  # convert from tau to beta
  # changed from fulltau to tau
  beta.foo <- transformUnconditional(parm = fulltau, 
    modelmatrix.int, data, from = "tau", to = "beta", 
    offset = offset)


  # set up for the bootstrap for the MLE
  theta.hat2 <- predict(model, model.type = "cond", 
    parm.type = "canon")
  theta.hat2 <- matrix(theta.hat2, nrow = n, ncol = nnode)
  est.tau <- matrix(nrow = npop, ncol = nboot) # changed from nrow = p
  b2 <- "try-error"; class(b2) <- "try-error"

  # set up for the bootstrap for the envelope
  theta.hat <- transformUnconditional(parm = fulltau, 
    modelmatrix.int, data, from = "tau", to = "theta", 
    offset = offset)
  theta.hat <- matrix(theta.hat, nrow = n, ncol = nnode)
  #theta.hat <- theta.hat2
  est  <- matrix(nrow = npop, ncol = nboot) # changed from nrow = p
  b <- "try-error"; class(b) <- "try-error"


  # construct matrix and array that specify
  # which components of varb correspond to 
  # Darwinian fitness
  amat.mat <- NULL
  if(class(amat) == "array"){
    amat.mat <- matrix(amat, nrow = npop, 
      byrow = TRUE)
  }
  if(class(amat) == "matrix"){
    amat.mat <- amat
    amat <- array(amat.mat, dim = c(npop, nnode, p))
  }


  # construct necessary quantities for the
  # hypothetical individuals in the newdata
  # argument
  if(is.null(newdata)){
    renewdata <- reshape(newdata, varying = list(vars),
      direction = "long", timevar = "varb",
      times = as.factor(vars), v.names = "resp")
    fit.renew <- as.numeric(grepl(fit.name, renewdata$varb))
    renewdata$fit <- fit.renew
    renewdata$root <- 1
  }

  modmat.renew <- modmat.new
  if(is.null(modmat.new)){
    modmat.renew <- model.matrix(formula, data = renewdata)
  }
  
  cond <- !(colnames(modmat.renew) %in% model$dropped)
  modmat.renew <- modmat.renew[, cond]
  M1.renew <- t(modmat.renew[,-index])
  M2.renew <- P %*% t(modmat.renew[,index])
  modmat.env.renew <- t(rbind(M1.renew,M2.renew))
  data.renew <- asterdata(newdata, vars = vars, pred = pred,
    group = rep(0, length(vars)), code = code, 
    families = families)
  origin.renew <- model$origin[1:npop,]
  offset.renew <- as.vector(origin.renew)


  # estimate expected Darwinian fitness using
  # the original data
  fit.env <- amat.mat %*% transformUnconditional(beta.foo,
    modmat.env.renew, data.renew, from = "beta", to = "mu", 
    offset = offset.renew)
  fit.MLE <- amat.mat %*% transformUnconditional(beta,
    modmat.renew, data.renew, from = "beta", to = "mu", 
    offset = offset.renew)


  # inital quantities to save computing time
  G.star <- matrix(0, nrow = k, ncol = u)
  P.star <- matrix(0, nrow = k, ncol = k)
  U.star <- U
  xstar <- xstar2 <- matrix(0, nrow = n, ncol = nnode)
  beta.env.star <- tau.env.star <- rep(0, p)
  target.star2 <- target.star <- rep(0, k)
  avar.star <- matrix(0, nrow = k, ncol = k)
  modelmatrix.star <- matrix(0, nrow = n*nnode, ncol = p)
  mu.star <- rep(0,nnode*n)
  aout4star <- aout4star2 <- model
  mu.renew <- mu.env.renew <- rep(0,npop)


  # the MLE bootstrap
  set.seed(13)
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
    #modelmatrix.star <- matrix(aout4star2$modmat, nrow = n*nnode)
    #mu.star <- predict(aout4star2, parm.type = "mean.value", 
    #  model.type = "unconditional")
    #tau.star <- crossprod(modelmatrix.star, mu.star)
    mu.renew <- transformUnconditional(aout4star2$coef, data.renew,
      modmat = modmat.renew, from = "beta", to = "mu",
      offset = offset.renew)
    est.tau[,iboot] <- amat.mat %*% mu.renew


    if(quiet == FALSE){
      if((iboot %% m) == 0){
        timer <- proc.time() - timer
        cat("iteration: ", iboot, " time: ", timer, "\n")
        timer <- proc.time()
      }
    }
  }


  # the envelope bootstrap
  set.seed(13)
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
    M <- t(matrix(modmat.model.int, nrow = n*nnode))
    avar.star <- (aout4star$fisher)[target.ind,target.ind]
    beta.env.star <- aout4star$coef

    if(method == "eigen"){
      G.star <- eigen(avar.star, symmetric = TRUE)$vec[,c(vectors)]  
      P.star <- tcrossprod(G.star)
    }
    if(method == "1d"){
      tau.env.star <- (M %*% as.vector(xstar))[target.ind]
      U.star <- tau.env.star %o% tau.env.star
      G.star <- manifold1Dplus(avar.star, U = U.star, u = u)
      P.star <- tcrossprod(G.star)
    }    
    M2 <- M[index,]; M2 <- P.star %*% M2
    M[index,] <- M2
    modelmatrix.star <- t(M)
    
    mu.star <- predict(aout4star, parm.type = "mean.value", 
      model.type = "unconditional")
    tau.env.star <- crossprod(modelmatrix.star, mu.star)
    beta.env.star <- transformUnconditional(parm = tau.env.star, 
      modelmatrix.star, data, from = "tau", to = "beta", 
      offset = offset)
    M1.renew <- t(modmat.renew[,-index])
    M2.renew <- P.star %*% t(modmat.renew[,index])
    modmat.env.renew <- t(rbind(M1.renew,M2.renew))      
    mu.env.renew <- transformUnconditional(parm = beta.env.star, 
      modmat.env.renew, data.renew, from = "beta", to = "mu", 
      offset = offset.renew)
    est[,iboot] <- amat.mat %*% mu.env.renew  


    #est[,iboot]  <- transformUnconditional(parm = aout4star$coef, 
    #  modelmatrix.star, data, from = "beta", to = "tau", 
    #  offset = offset, tolerance = 1e-200)
    

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
  table <- cbind(fit.env, sqrt(diag(S)), fit.MLE, 
    sqrt(diag(S2)), ratio)
  colnames(table) <- c("env","se(env)","MLE","se(MLE)","ratio")  
  out <- list(u = u, table = table, S = S, S2 = S2, 
    env.boot.out = est, MLE.boot.out = est.tau)
  return(out)
  
}