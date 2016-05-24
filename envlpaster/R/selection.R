

selection <- function(parm, index, model, data = NULL, 
  alpha = 0.05, type = c("canonical","mean-value"), 
  method = c("eigen","1d")){

  # important initial quantities
  modmat <- model$modmat
  offset <- as.vector(model$origin)
  nind <- dim(modmat)[1]
  nnode <- dim(modmat)[2]
  m <- length(index)
  p <- length(parm)
  modmat.mat <- matrix(modmat, nrow = nind * nnode)
  x <- model$x
  root <- model$root
  pred <- model$pred
  fam <- model$fam

  # get the spectral structure of FI or inv FI
  avar <- model$fisher
  if(type == "canonical") avar <- solve(avar, tol = 1e-20)
  eig.avar.targ <- eigen(avar[index,index])

  # get important quantities for 1D algorithm
  avar.targ <- avar[index, index]
  U <- parm %*% t(parm)
  U.targ <- U[index, index]
  
  # initialize internal quantities
  beta.cand.test <- NULL
  beta <- model$coef
  beta.foo <- rep(0, p)
  tau.int <- parm
  aic.env <- bic.env <- pval <- 0
  
  # consider every possible envelope estimator constructed 
  # using eigenspaces
  if(method == "eigen"){
    beta.cand.test <- matrix(unlist(lapply(1:m, FUN = function(j) {
      cand <- combs(1:m, j)
      matrix(unlist(lapply(1:nrow(cand), FUN = function(k) {

        # obtain the eigenspace of interest and the projection
        # into that eigenspace
        G <- eig.avar.targ$vec[,cand[k,]]
        P <-  projection(G)
      
        # build the model matrix Menv
        M <- t(modmat.mat)
        #M2 <- M[(1:p %in% index),]; M2 <- P %*% M2
        M[(1:p %in% index),] <- P %*% M[(1:p %in% index),]
        modmat.mat.env <- t(M)
        modmat.int <- array(modmat.mat.env, c(nind, nnode, p))

        # obtain the corresponing beta vector
        if(type == "mean-value"){
          tau.env.int <- crossprod(P, parm[index]) 
          tau.int[index] <- tau.env.int
          ind.int <- max(which(!1:p %in% index)) + length(cand[k,])
          beta.foo[1:(ind.int)] <- suppressWarnings(try(
            aster(x, root, pred, fam, modmat = modmat.int, 
              origin = model$origin, maxiter = 10000)$coef, 
          silent = TRUE))
          #beta.foo <- try(transformUnconditional(parm = tau.int, 
          #  modmat.mat.env, data, from = "tau", to = "beta", 
          #  offset = offset), silent = TRUE)
          #print(beta.foo)
        }
        
        if(type == "canonical"){
          beta.foo[index] <- P %*% beta.foo[index]
        }

        bic.env <- aic.env <- 1e9
        pval <- 0
        if(class(beta.foo) != "try-error"){
          # obtain selection criteria values
          env <-mlogl(beta.foo, pred, fam, x, root, 
            modmat = modmat.int, type = "unconditional")$value
          full <-mlogl(beta, pred, fam, x, root, 
            modmat = modmat, type = "unconditional")$value
          df <- j*(m-j) + j + j*(j+1)/2 + (m-j)*(m-j+1)/2
          df.full <- m + m*(m+1)/2
          LRT <- 2*(env - full)
          if(j < m) pval <- pchisq(LRT, df = df.full - df, lower.tail = F)
          if(j == m) pval <- 1
          bic.env <-2*env + df*log(nind)
          aic.env <- 2*df + 2*env
        }
        
        things <- c(j, aic.env, bic.env, pval)
        return(things)
      } )))
    } )), ncol = 4, byrow = TRUE)
  }

  # consider every possible envelope estimator using 
  # the 1D algorithm
  if(method == "1d"){
    beta.cand.test <- matrix(unlist(lapply(1:m, FUN = function(j) {
      if(j < m){

        # obtain the eigenspace of interest and the projection
        # into that eigenspace      
        G <- manifold1Dplus(M = avar.targ, U = U.targ, u = j)
        P <- projection(G)

        # build the model matrix Menv
        M <- t(modmat.mat)
        M2 <- M[(1:p %in% index),]; M2 <- P %*% M2
        M[(1:p %in% index),] <- M2
        modmat.mat.env <- t(M)
        modmat.int <- array(modmat.mat.env, c(nind, nnode, p))        

        # obtain the corresponing beta vector
        if(type == "mean-value"){
          tau.env.int <- crossprod(P, parm[index]) 
          tau.int[index] <- tau.env.int
          beta.foo <- transformUnconditional(parm = tau.int, 
            modmat.mat.env, data, from = "tau", to = "beta", 
            offset = offset, tolerance = 1e-200)
        }

        if(type == "canonical"){
          beta.foo[index] <- P %*% beta.foo[index]
        }

        # obtain selection criteria values  
        env <-mlogl(beta.foo, pred, fam, x, root, 
          modmat = modmat.int, type = "unconditional")$value
        full <-mlogl(beta, pred, fam, x, root, 
          modmat = modmat, type = "unconditional")$value
        k <- j*(m-j) + j + j*(j+1)/2 + (m-j)*(m-j+1)/2
        k.full <- m + m*(m+1)/2
        LRT <- 2*(env - full)
        pval <- pchisq(LRT, df = k.full - k, lower.tail = F)
        bic.env <-2*env + k*log(nind)
        aic.env <- 2*k + 2*env
      }
  
      if(j == m){
        full <- mlogl(beta, pred, fam, x, root, 
          modmat = modmat, type = "unconditional")$value
        k <- (m + m*(m+1)/2)
        bic.env <-2*full + k*log(nind) 
        aic.env <- 2*k + 2*full
        pval <- 1
      }  
  
      things <- c(j, aic.env, bic.env, pval)
      return(things)
    })), ncol = 4, byrow = TRUE)
  }


  # select the dimension for each of the criteria
  aic <- bic <- LRT <- m
  if(method == "1d"){
    aic <- which(beta.cand.test[,2] == min(beta.cand.test[,2]))
    bic <- which(beta.cand.test[,3] == min(beta.cand.test[,3]))
    LRT <- min(which(beta.cand.test[,4] > alpha))
  }

  # select the eigenspace for each of the criteria
  if(method == "eigen"){
    dims <- nrow(combs(1:m,1))   
    for(u in 2:m) dims[u] <- nrow(combs(1:m, u)) + dims[u-1]     
    aic <- which(beta.cand.test[,2] == min(beta.cand.test[,2]))
    bic <- which(beta.cand.test[,3] == min(beta.cand.test[,3]))
    LRT <- min(which(beta.cand.test[,4] > alpha))
    u.aic <- min(which(dims >= aic))
    u.bic <- min(which(dims >= bic))
    u.LRT <- min(which(dims >= LRT))
    if(u.aic == 1) aic <- combs(1:m, u.aic)[c(aic), ]    
    if(u.aic >= 2) aic <- combs(1:m, u.aic)[c(aic - dims[u.aic - 1]), ]
    if(u.bic == 1) bic <- combs(1:m, u.bic)[c(bic), ]    
    if(u.bic >= 2) bic <- combs(1:m, u.bic)[c(bic - dims[u.bic - 1]), ]
    if(u.LRT == 1) LRT <- combs(1:m, u.LRT)[c(LRT), ]    
    if(u.LRT >= 2) LRT <- combs(1:m, u.LRT)[c(LRT - dims[u.LRT - 1]), ]
  }

  beta.cand.test <- as.data.frame(beta.cand.test)
  colnames(beta.cand.test) <- c("u","aic","bic","p-val")

  out <- list(aic = aic, bic = bic, LRT = LRT, out = beta.cand.test)
  return(out)
}


