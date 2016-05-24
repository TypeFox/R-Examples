



secondboot <- function(k, nboot2, out, model, index, data, amat, 
    newdata, method = c("eigen","1d"))
{

  # extract necessary components from the top-level of
  # bootstrapping 
  env.boot.out <- out$env.boot.out
  env.1d.boot.out <- out$env.1d.boot.out
  MLE.boot.out <- out$MLE.boot.out
  MLE.tau.boot <- out$MLE.tau.boot
  env.tau.boot <- out$env.tau.boot
  env.1d.tau.boot <- out$env.1d.tau.boot
  P.list <- out$P.list
  P.1d.list <- out$P.1d.list
  vectors.list <- out$vectors.list
  u.1d.list <- out$u.1d.list
  nboot <- ncol(env.boot.out)
  npop <- nrow(env.boot.out)


  # specify necessary quantities for secondboot function
  # not extracted in the above
  aout4star2 <- b2 <- model
  beta <- model$coef
  p <- length(beta)
  n <- nrow(model$x)
  nnode <- ncol(model$x)
  modmat.mat <- matrix(model$modmat, nrow = n * nnode)
  mu <- predict(model, parm.type = "mean.value", 
    model.type = "unconditional")
  #tau <- crossprod(modmat.mat, mu)
  tau <- MLE.tau.boot[, k]
  offset <- as.vector(model$origin)
  code <- data$code
  families <- data$families
  vars <- colnames(model$x)
  fam <- model$fam
  pred <- model$pred
  root <- model$root  


  # necessary for trial run
  cov <- V <- var.Efron <- NULL
  MLE.tau.boot.subsample <- NULL
  est.env.subsample <- NULL


  # the array (matrix) that specifies Darwinian fitness
  amat.mat <- NULL
  if(class(amat) == "array"){
    amat.mat <- matrix(amat, nrow = npop, 
    byrow = TRUE)
  }
  if(class(amat) == "matrix"){
    amat.mat <- amat
    amat <- array(amat.mat, dim = c(npop, nnode, p))
  }

  # initialize important quantities
  P.foo <- M.foo <- fit.foo <- foo <- NULL
  vectors.foo <- u <- NULL

  if(method == "eigen"){
    foo <- env.tau.boot[, k] 
    vectors.foo <- vectors.list[[k]]
    u <- length(vectors.foo)
    P.foo <- P.list[[k]] 
    M.foo <- modmat.mat; M.foo[, index] <- M.foo[, index] %*% P.foo
    fit.foo <- rowSums(out$env.boot.out) / nboot
  }

  if(method == "1d"){
    u <- u.1d.list[[k]]
    fit.foo <- rowSums(out$env.1d.boot.out) / nboot
    if(u == length(index)){
      foo <- MLE.tau.boot[, k]
      P.foo <- P.1d.list[[k]] 
      M.foo <- modmat.mat
    }
    if(u < length(index)){
      foo <- env.1d.tau.boot[, k] 
      P.foo <- P.1d.list[[k]] 
      M.foo <- modmat.mat; M.foo[, index] <- M.foo[, index] %*% P.foo
    }
  }

  # theta values
  theta.samp <- transformUnconditional(foo, M.foo, data, 
    from = "tau", to = "theta", offset = offset)
  theta.samp <- matrix(theta.samp, nrow = n, ncol = nnode)
  
  # setup of model matrices and asterdata object 
  # obtained from top level of bootstrapping
  cond <- !(colnames(modmat.renew) %in% model$dropped)
  modmat.renew <- modmat.env.renew <- modmat.renew[, cond]
  modmat.env.renew[, index] <- 
    modmat.env.renew[, index] %*% P.foo
  data.renew <- asterdata(newdata, vars = vars, pred = pred,
    group = rep(0, length(vars)), code = code, 
    families = families)
  origin.renew <- model$origin[1:npop,]
  offset.renew <- as.vector(origin.renew)    

  # initial quantities
  #M.foo.array <- array(M.foo, dim = dim(model$modmat))
  #print(dim(M.foo.array))
  #print(dim(model$modmat))
  r <- length(index)
  U.2nd <- P.2nd <- Sigma.uu <- matrix(0, r, r)
  Gamma.2nd <- matrix(0, nrow = r, ncol = u)
  M.2nd <- M.foo
  est.env.subsample <- matrix(0, nrow = npop, ncol = nboot2)
  MLE.tau.boot.subsample <- matrix(0, nrow = p, ncol = nboot2)
  sd.Efron <- 0

  # generate many new resamples    
  for(j in 1:nboot2){

    xstar.samp <- raster(theta.samp, pred, fam, root)
    
    # fit the aster model using the secondary generated data
    class(b2) <- class(try(aout4star2 <- aster(xstar.samp, root, 
      pred, fam, model$modmat, parm = beta), silent = TRUE))[1]

    if(class(b2) == "try-error"){
      class(b2) <- class(try(aout4star2 <- aster(xstar.samp, root, pred,
        fam, model$modmat, parm = beta, 
        method = "nlm"), silent = TRUE))[1] 
    }

    if(class(b2) == "try-error"){
      class(b2) <- class(try( aout4star2 <- aster(xstar.samp, root, 
        pred, fam, model$modmat), silent = TRUE))[1]
    }

    if(class(b2) == "try-error"){
      class(b2) <- class(try(aout4star2 <- aster(xstar.samp, root, pred,
        fam, model$modmat, method = "nlm"), silent = TRUE))[1] 
    }


    # get tau and mu from this fit
    mu.star <- predict(aout4star2, parm.type = "mean.value", 
      model.type = "unconditional")

    #tau.renew <- transformUnconditional(aout4star2$coef, data,
    #  modmat = modmat.mat, from = "beta", to = "tau",
    #  offset = offset)      

    # get beta and tau using the model matrix containing
    # the projection into the envelope
    Sigma.uu <- aout4star2$fisher[index, index]
    if(method == "eigen"){
      if(length(vectors.foo) < length(index)){
        Gamma.2nd <- eigen(Sigma.uu, symmetric = TRUE)$vec[, vectors.foo]
        P.2nd <- projection(Gamma.2nd)
        M.2nd[, index] <- M.2nd[, index] %*% P.2nd
      }
      if(length(vectors.foo) == length(index)){
        M.2nd <- M.foo
      }
    }
    if(method == "1d"){
      if(u < length(index)){
        tau.star <- crossprod(M.foo, mu.star)
        U.2nd <- tau.star[index] %*% t(tau.star[index]) 
        Gamma.2nd <- manifold1Dplus(Sigma.uu, U = U.2nd, u = u)
        P.2nd <- tcrossprod(Gamma.2nd)
      }
      if(u == length(index)){ 
        P.2nd <- diag(length(index))
        M.2nd[, index] <- M.foo[, index] %*% P.2nd  
      }
    }


    tau.env.star <- crossprod(M.2nd, mu.star)
    beta.env.star <- transformUnconditional(parm = tau.env.star, 
      M.2nd, data, from = "tau", to = "beta", 
      offset = offset)


    # compute the envelope estimator of expected 
    # Darwinian fitness
    modmat.env.renew[, index] <- modmat.renew[,index] %*% P.2nd
    mu.env.renew <- transformUnconditional(parm = beta.env.star, 
      modmat.env.renew, data.renew, from = "beta", to = "mu", 
      offset = offset.renew)
    est.env.subsample[, j] <- (amat.mat %*% mu.env.renew)
    M.2nd <- M.foo

    # store the canonical statistic value
    MLE.tau.boot.subsample[, j] <- tau.env.star

    # compute the Efron sd estimator99i5
    if(j == nboot2){ 

      B <- t(MLE.tau.boot.subsample - rowSums(MLE.tau.boot.subsample)/nboot2)
      V <- (t(B) %*% B) / nboot2;  eig.V <- eigen(V)
      env.centered <- t(est.env.subsample - fit.foo)
      cov <- t(B) %*% env.centered / nboot2

      var.Efron <- t(cov) %*% eig.V$vec %*% diag(1/eig.V$val) %*% 
        t(eig.V$vec) %*% cov
      sd.Efron <- sqrt(diag(var.Efron))
    }
  }

  ### output 
  out <- list(sd.Efron = sd.Efron, cov = cov, V = V, 
    MLE.tau.boot.subsample = MLE.tau.boot.subsample, 
    est.env.subsample = est.env.subsample)
  return(out)
}





