

targetbootcanon <- function(model, nboot, index, u, 
  code, families, quiet = FALSE, m = 100){
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
  dimensions <- dim(modmat.model)
  fam <- model$fam
  pred <- model$pred
  root <- model$root
  x <- model$x
  vars <- colnames(x)  
  n <- nrow(x)
  nnode <- ncol(x)
  modmat.mat <- matrix(modmat.model, nrow = n * nnode)
  origin <- model$origin
  offset <- as.vector(origin)

  # run the 1-d algorithm w.r.t. the target
  # construct the envelope estimator
  avar <- solve(model$fisher, symmetric = TRUE)[target.ind,target.ind]
  foo <- manifold1Dplus(M = avar, U = U, u = u)
  P <- projection(foo)
  beta.env <- P %*% target  
  fullbeta <- rep(0,p)
  fullbeta[target.ind] <- beta.env
  fullbeta[nuis.ind] <- nuisance

  # change the model matrix
  M <- t(modmat.mat)
  M2 <- M[index,]; M2 <- P %*% M2
  M[index,] <- M2
  modmat.mat.env <- t(M)
  modmat.mat.int <- array(modmat.mat.env, dimensions)


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
  P.star <- fisher.star <- foot.star <- 
    matrix(0, nrow = p, ncol = p)


  # the MLE bootstrap
  for(iboot in 1:nboot){

    xstar2 <- raster(theta.hat2, pred, fam, root)
    colnames(xstar2) <- vars
    
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

    # fix the data reference, we may need a data
    # frame argument in order to procede
    xstar <- raster(theta.hat, pred, fam, root)
    colnames(xstar) <- vars
    vec <- as.vector(xstar)
    tau <- crossprod(modmat.mat.env, vec)
    data.star <- data
    print(vars)
    print(colnames(data))
    cond <- colnames(data) %in% vars
    data.star[,cond] <- xstar
    aster.data <- asterdata(data.star, vars, pred, 
      group = rep(0,nnode), code = code, delta = rep(0,nnode),
      families = families)

    beta.int <- transformUnconditional(tau, modmat.mat.env,
      aster.data, from = "tau", to = "beta", offset = offset)
    jacob.env <- jacobian(beta.int, aster.data, from = "beta", 
      to = "tau", transform = "unconditional", 
      modmat = modmat.mat, offset = offset)
    foot.star[nuis.ind,nuis.ind] <- diag(length(nuis.ind))
    foot.star[target.ind,target.ind] <- P
    avar.star <- foot.star %*% jacob.env %*% foot.star
    G.star <- eigen(avar.star[target.ind,target.ind])$vec[,1]
    P.star <- tcrossprod(G.star)
    est[,iboot] <- P.star %*% beta.int[target.ind]


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
  S <- var(t(est))
  S2 <- var(t(est.beta))
  ratio <- sqrt(diag(S2) / diag(S))
  table <- cbind(beta.env, sqrt(diag(S)), beta[target.ind], 
    sqrt(diag(S2)), ratio)
  colnames(table) <- c("env","se(env)","MLE","se(MLE)","ratio")
  out <- list(u = u, table = table, S = S, S2 = S2, env.boot.out = est,
    MLE.boot.out = est.beta)
  return(out)
}
##########################################################
















eigenbootcanon <- function(model, nboot, index, vectors, 
  code, families, quiet = FALSE, m = 100){
  timer <- proc.time()

  # extract important envelope initial quantities
  # obtain the target and nuisance parameters
  # build U w.r.t. the target parameters
  u <- length(vectors)
  beta <- model$coef
  p <- length(beta)
  nuis.ind <- c(1:p)[-c(index)]
  target.ind <- index
  k <- length(index)
  nuisance <- beta[nuis.ind]


  # important aster model quantities
  x <- model$x
  vars <- colnames(x)  
  n <- nrow(x)
  nnode <- ncol(x)
  formula <- model$formula
  modmat.model <- model$modmat
  dimensions <- dim(modmat.model)
  modmat.mat <- matrix(modmat.model, nrow = n * nnode)
  fam <- model$fam
  pred <- model$pred
  root <- model$root
  offset <- as.vector(model$origin)


  # run the eigenspace approach w.r.t. the target
  # construct the envelope estimator
  avar <- solve(model$fisher, symmetric = TRUE)[target.ind,target.ind]
  eig <- eigen(avar, symmetric = TRUE)  
  G <- eig$vec[,c(vectors)]
  P <- tcrossprod(G)
  beta.env <- crossprod(P,beta[target.ind])
  fullbeta <- rep(0,p)
  fullbeta[target.ind] <- beta.env
  fullbeta[nuis.ind] <- nuisance


  # change the model matrix
  M <- t(modmat.mat)
  M2 <- M[index,]; M2 <- P %*% M2
  M[index,] <- M2
  modmat.mat.env <- t(M)
  modmat.mat.int <- array(modmat.mat.env, dimensions)


  # set up for the bootstrap for the envelope
  model2 <- model; class(model2) <- c("aster")
  theta.hat <- predict(model2, x, model$root, 
    model.type = "cond", parm.type = "canon",
    newcoef = fullbeta, modmat = modmat.mat.int)
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
  aster.data <- asterdata(data, vars, pred, 
    group = rep(0,nnode), code = code, delta = rep(0,nnode),
    families = list("bernoulli", "zero.truncated.poisson"))
  data.star <- data
  xstar <- xstar2 <- matrix(0, nrow = n, ncol = nnode)
  target.star2 <- beta.env.star <- target.star <- rep(0, k)
  avar.star <- matrix(0, nrow = k, ncol = k)
  G.star <- matrix(0, nrow = p, ncol = u)
  P.star <- fisher.star <- foot.star <- 
    matrix(0, nrow = p, ncol = p)


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
    vec <- as.vector(xstar)
    tau <- crossprod(modmat.mat.env, vec)
    data.star <- data
    cond <- colnames(data) %in% vars
    data.star[,cond] <- xstar
    aster.data <- asterdata(data.star, vars, pred, 
      group = rep(0,nnode), code = code, delta = rep(0,nnode),
      families = families)

    beta.int <- transformUnconditional(tau, modmat.mat.env,
      aster.data, from = "tau", to = "beta", offset = offset)
    jacob.env <- jacobian(beta.int, aster.data, from = "beta", 
      to = "tau", transform = "unconditional", 
      modmat = modmat.mat, offset = offset)
    foot.star[nuis.ind,nuis.ind] <- diag(length(nuis.ind))
    foot.star[target.ind,target.ind] <- P
    avar.star <- foot.star %*% jacob.env %*% foot.star
    G.star <- eigen(avar.star[target.ind,target.ind])$vec[,1]
    P.star <- tcrossprod(G.star)
    est[,iboot] <- P.star %*% beta.int[target.ind]

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
  S <- var(t(est))
  S2 <- var(t(est.beta))
  ratio <- sqrt(diag(S2) / diag(S))
  table <- cbind(beta.env, sqrt(diag(S)), beta[target.ind], 
    sqrt(diag(S2)), ratio)
  colnames(table) <- c("env","se(env)","MLE","se(MLE)","ratio")
  out <- list(u = u, table = table, S = S, S2 = S2, env.boot.out = est,
    MLE.boot.out = est.beta)    
  return(out)
}
##########################################################

