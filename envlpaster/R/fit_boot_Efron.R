


fit.boot.Efron <- function(model, nboot, index, vectors = NULL, dim = NULL,
    data, amat, newdata, modmat.new = NULL, renewdata = NULL, 
    criterion = c("AIC","BIC","LRT"), alpha = 0.05, fit.name = NULL, 
    method = c("eigen","1d"), quiet = FALSE)
{

  # stopping condition
  #stopifnot(!(is.null(vectors) & is.null(u)))
  #stopifnot(!(is.null(fit.name) & is.null(renewdata)))
  #timer <- proc.time()

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
  modmat.mat <- matrix(modmat.model, nrow = n * nnode)
  dimensions <- dim(modmat.model)
  offset <- as.vector(model$origin)
  
  # obtain tau
  mu <- predict(model, parm.type = "mean.value", 
    model.type = "unconditional")
  tau <- crossprod(modmat.mat, mu)
  nuis.ind <- c(1:p)[!index]
  k <- length(index)
  target <- tau[index]

  # obtain the projection used to construct the envelope 
  # estimator
  avar <- (model$fisher)[index,index]
  U <- P <- NULL; u <- length(index)
  if(method == "eigen"){
    u <- length(vectors)
    eig <- eigen(avar, symmetric = TRUE)  
    G <- eig$vec[,c(vectors)]
    P <- tcrossprod(G)  
  }
  if(method == "1d"){
    u <- dim
    U <- target %o% target
    G <- manifold1Dplus(M = avar, U = U, u = u)
    P <- tcrossprod(G)
  }

  tau.env <- crossprod(P,target)
  fulltau <- tau
  fulltau[index] <- tau.env

  # change the model matrix for the envelope estimator
  modelmatrix.int <- modmat.mat
  modelmatrix.int[, index] <- modelmatrix.int[, index] %*% P
  modmat.model.int <- array(modmat.mat, dimensions)

  # convert from tau to beta changed from fulltau to tau
  beta.foo <- transformUnconditional(parm = fulltau, 
    modelmatrix.int, data, from = "tau", to = "beta", 
    offset = offset)

  # set up for the bootstrap for the MLE
  theta.hat2 <- predict(model, model.type = "cond", 
    parm.type = "canon")
  theta.hat2 <- matrix(theta.hat2, nrow = n, ncol = nnode)
  MLE.tau.boot <- matrix(nrow = npop, ncol = nboot) # changed from nrow = p
  b2 <- "try-error"; class(b2) <- "try-error"

  # set up for the bootstrap for the envelope estimator
  # (this envelope estimator corresponds to the selected method)
  theta.hat <- transformUnconditional(parm = fulltau, 
    modelmatrix.int, data, from = "tau", to = "theta", 
    offset = offset)
  theta.hat <- matrix(theta.hat, nrow = n, ncol = nnode)
  #theta.hat <- theta.hat2
  est  <- matrix(nrow = npop, ncol = nboot) # changed from nrow = p
  b <- "try-error"; class(b) <- "try-error"


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


  # initial quantities
  est <- est2 <- est.1d <- matrix(0, nrow = npop, ncol = nboot)
  MLE.tau.boot <- env.tau.boot <- env.1d.tau.boot <- matrix(0, nrow = p, ncol = nboot)
  P.list <- P.1d.list <-  NULL
  class(P.list) <- class(P.1d.list) <- "list"
  vectors.list <- u.1d.list <- NULL
  class(vectors.list) <- class(u.1d.list) <- "list"
  out <- NULL
  table <- NULL


  #modmat.renew <- model.matrix(m1$formula, data = renewdata)
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



  # the parametric bootstrap which takes Efron's procedure into account
  for(k in 1:nboot){

    # generate a resample of responses from the MLE
    xstar2 <- raster(theta.hat2, pred, fam, root)

    # fit the aster model to the resampled data obtained from 
    # the MLE of theta
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


    # MLE of expected Darwinian fitness
    tau.renew <- transformUnconditional(aout4star2$coef, data,
      modmat = modmat.mat, from = "beta", to = "tau",
      offset = offset)
    phi.MLE.renew <- offset.renew + modmat.renew %*% aout4star2$coef
    mu.MLE.renew <- transformSaturated(parm = phi.MLE.renew, 
      data.renew, from = "phi", to = "mu")
    MLE.tau.boot[, k] <- tau.renew
    Dar.fit.MLE <- amat.mat %*% mu.MLE.renew


    # generate a resample of responses from the 
    # envelope estimator
    xstar <- raster(theta.hat, pred, fam, root)

    # fit the aster model to the resampled data obtained from 
    # the envelope estimator of theta
    class(b2) <- class(try(aout4star2 <- aster(xstar, root, pred,
      fam, modmat.model.int, parm = beta), silent = TRUE))[1] 
    
    if(class(b2) == "try-error"){
      class(b2) <- class(try(aout4star2 <- aster(xstar, root, pred,
        fam, modmat.model.int, parm = beta, 
        method = "nlm"), silent = TRUE))[1] 
    }

    if(class(b2) == "try-error"){
      class(b2) <- class(try( aout4star2 <- aster(xstar, root, 
        pred, fam, modmat.model.int), silent = TRUE))[1]
    }

    if(class(b2) == "try-error"){
      class(b2) <- class(try(aout4star2 <- aster(xstar, root, pred,
        fam, modmat.model.int, method = "nlm"), silent = TRUE))[1] 
    }


    # MLE of expected Darwinian fitness
    tau.renew <- transformUnconditional(aout4star2$coef, data,
      modmat = modelmatrix.int, from = "beta", to = "tau",
      offset = offset)
    #phi.MLE.renew <- offset.renew + modmat.renew %*% aout4star2$coef
    #mu.MLE.renew <- transformSaturated(parm = phi.MLE.renew, 
    #  data.renew, from = "phi", to = "mu")

    # selection of the envelope dimension using the 1d algorithm
    barbaz <- NULL
    dim.1d <- 0
    if(method == "1d"){
      barbaz <- selection(tau.renew, index, model = aout4star2, data = data, 
        alpha = alpha, type = "mean-value", method = "1d")
      if(criterion == "AIC") dim.1d <- barbaz$aic
      if(criterion == "BIC") dim.1d <- barbaz$bic
      if(criterion == "LRT") dim.1d <- barbaz$LRT
    }

    u.1d.list[[k]] <- dim.1d

    # selection of the envelope eigenstructure
    fubar <- NULL
    vectors <- NULL
    if(method == "eigen"){
      fubar <- selection(tau.renew, index, model = aout4star2, data = data, 
        alpha = alpha, type = "mean-value", method = "eigen")
      if(criterion == "AIC") vectors <- fubar$aic
      if(criterion == "BIC") vectors <- fubar$bic
      if(criterion == "LRT") vectors <- fubar$LRT    
    }

    vectors.list[[k]] <- vectors

    if(!quiet){
      cat("iteration: " , k, " indices: ", vectors, " dim.1d: ", dim.1d, "\n")
    }

    # get the bootstrapped envelope estimators    
    M <- matrix(modmat.model.int, nrow = n * nnode)
    avar.star <- (aout4star2$fisher)[index,index]
    beta.env.star <- aout4star2$coef
    mu.star <- predict(aout4star2, parm.type = "mean.value", 
      model.type = "unconditional")
    
    # envelope estimator of expected Darwinian fitness using
    # eigenstructures
    Dar.fit.env <- Dar.fit.MLE
    if(method == "eigen"){    
      G.star <- eigen(avar.star, symmetric = TRUE)$vec[,c(vectors)]  
      P.star <- tcrossprod(G.star)
      P.list[[k]] <- P.star
      modelmatrix.star <- M
      modelmatrix.star[, index] <- modelmatrix.star[, index] %*% P.star

      tau.env.star <- crossprod(modelmatrix.star, mu.star)
      env.tau.boot[, k] <- tau.env.star
      beta.env.star <- transformUnconditional(parm = tau.env.star, 
        modelmatrix.star, data, from = "tau", to = "beta", 
        offset = offset)
      M1.renew <- t(modmat.renew[,-index])
      M2.renew <- P.star %*% t(modmat.renew[,index])
      modmat.env.renew <- t(rbind(M1.renew,M2.renew))      
      mu.env.renew <- transformUnconditional(parm = beta.env.star, 
        modmat.env.renew, data.renew, from = "beta", to = "mu", 
        offset = offset.renew)
      Dar.fit.env <- amat.mat %*% mu.env.renew  
    }

    # envelope estimator of expected Darwinian fitness using the
    # 1d algorithm
    Dar.fit.env.1d <- Dar.fit.MLE  
    if(method == "1d"){
      if(dim.1d < length(index)){
        up.env.star <- crossprod(M, mu.star)[index]
        U.1d.star <- up.env.star %o% up.env.star
        G.1d.star <- manifold1Dplus(avar.star, U = U.1d.star, u = dim.1d)
        P.1d.star <- tcrossprod(G.1d.star)
        P.1d.list[[k]] <- P.1d.star

        modelmatrix.star <- M
        modelmatrix.star[, index] <- modelmatrix.star[, index] %*% P.1d.star

        tau.env.star <- crossprod(modelmatrix.star, mu.star)
        env.1d.tau.boot[, k] <- tau.env.star
        beta.env.star <- transformUnconditional(parm = tau.env.star, 
          modelmatrix.star, data, from = "tau", to = "beta", 
          offset = offset)
        modmat.env.renew <- modmat.renew
        modmat.env.renew[, index] <- modmat.env.renew[, index] %*% P.1d.star
        mu.env.renew <- transformUnconditional(parm = beta.env.star, 
          modmat.env.renew, data.renew, from = "beta", to = "mu", 
          offset = offset.renew)
        Dar.fit.env.1d <- amat.mat %*% mu.env.renew 
      }

      if(dim.1d == length(index)){ 
        P.1d.list[[k]] <- diag(length(index))
        env.1d.tau.boot[, k] <- tau.renew
      }
    }


    # the stored expected Darwinian fitness estimates
    est[, k] <- Dar.fit.env
    est2[, k] <- Dar.fit.MLE
    est.1d[, k] <- Dar.fit.env.1d


    if(k == nboot){
      #means <- apply(est, FUN = mean, MARGIN = 1)
      #S <- var(t(est)); S2 <- var(t(est2)); S.1d <- var(t(est.1d))
      #ratio <- sqrt( diag(S2) / diag(S) )
      #table <- cbind(Dar.fit.env, sqrt(diag(S)), Dar.fit.MLE, 
      #  sqrt(diag(S2)), ratio)
      #colnames(table) <- c("env","se(env)","MLE","se(MLE)","ratio")  

      #ratio.1d <- sqrt( diag(S.1d) / diag(S) )
      #table.1d <- cbind(Dar.fit.env, sqrt(diag(S)), Dar.fit.env.1d, 
      #  sqrt(diag(S.1d)), ratio.1d)
      #colnames(table.1d) <- c("env","se(env)","env.1d","se(env.1d)","ratio")  
      out <- list( env.boot.out = est, MLE.boot.out = est2, 
        env.1d.boot.out = est.1d, MLE.tau.boot = MLE.tau.boot, 
        env.tau.boot = env.tau.boot, env.1d.tau.boot = env.1d.tau.boot, 
        P.list = P.list, P.1d.list = P.1d.list,
        vectors.list = vectors.list, u.1d.list = u.1d.list)
    }
  }

  return(out)
}





