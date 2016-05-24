# Construct a more efficient loop that considers the most important
# eigenvectors in that order
scanner <- function(M, coef, u){
  
  # initial values
  p <- ncol(M)
  indicesk <- par <- final <- values <- indices <- candidate <- NULL
  H <- as.matrix(rep(0,p))
  K <- as.matrix(rep(0,p))
  min <- 0
  out <- list()
  eigen.out <- eigen(M, symmetric = TRUE)
  eigenvect <- eigen.out$vectors
  eigenvalues <- eigen.out$values
  
  # find the most influential eigenvector
  internal <- rep(0,p)
  for(j in 1:p){

    # the objective function to minimize
    H <- eigenvect[,c(j)]
    f <- function(w){
      internal <- w * H
      y <- (coef - internal)
      x <- (t(y) %*% y)^5
      x
    }

    init <- 0
    error <- try( foo <- optim(init, fn = f, method = "Brent", 
	  lower = -1000000, upper = 1000000) )
    if(class(error) != "try-error"){    
      candidate[j] <- foo$value
      values[j] <- foo$par
    }
  }

  indices[1] <- min <- which(min(candidate) == candidate)
  final[1] <- values[min]
  
  if(u == 1){
    vec <- as.matrix(eigenvect[,c(indices)] * final[1])
  }


  # based on the selection of the most influential eigenvector,
  # find the remaining eigenvectors
  if(u > 1){
    for(j in 2:u){
    
      # inital values
      candidate <- values <- NULL 
      H <- as.matrix(eigenvect[,c(indices)])
      colnames(H) <- as.character(c(indices))
      K <- as.matrix(eigenvect[,-c(indices)])
      indicesk <- as.numeric(subset(t(as.matrix(c(1:p))), subset = TRUE, select = -c(indices)))
      colnames(K) <- as.character(c(indicesk))
      lenk <- ncol(K)
      lenh <- ncol(H)

      # an internal loop to pick the next eigenvector
      for(k in 1:lenk){
        f <- function(w){
          for(h in 1:lenh){
            if(h == 1) internal <- final[1] * H[,c(1)]
            else if(h > 1) internal <- internal + final[h] * H[,c(h)]
          }
          internal <- internal + w * K[,c(k)]
          y <- (coef - internal)
          x <- (t(y) %*% y)^5
          x
        }

        init <- 0
        error <- try( foo <- optim(init, fn = f, method = "Brent", 
	      lower = -1000000, upper = 1000000) )  
        if(class(error) != "try-error"){    
           candidate[k] <- foo$value
           values[k] <- foo$par
        }
      }
      min <- which(min(candidate) == candidate)
      indices[j] <- indicesk[min]
      final[j] <- values[min]
    }

    H <- as.matrix(eigenvect[,c(indices)])
    colnames(H) <- as.character(c(indices))

    vec <- rep(0,p)
    for(j in 1:u) vec <- vec + final[j] * H[,j]

  }

  table <- cbind(vec, coef)
  colnames(table) <- c("estimated", "provided")
  prop <- sum(eigenvalues[-c(indices)]) / sum(eigenvalues)

  out = list(indices = indices, table = table, G = H,
    prop = prop)
  out  

}




##########################################################
# diagnostics (canonical parameterization)
#
# arguments:
# model - an aster model fit
# modelmat - model matrix corresponding to the aster model
# u - proposed dimension of the envelope
#
# output:
# table - an output table containing the aster model MLE for
#         beta, the envelope estimator, which components of
#         the envelope estimator are w/n two SEs of the MLE,
#         and the efficiency gains presented as a ratio 
# pval - p value from the LRT (asymptotic reference distn
#        is assumed) comparing the envelope model (H0) to
#        the aster model (Ha)
# aic - returns the selected model based on aic values
# bic - returns the selected model based on bic values
##########################################################
#diagnostics <- function(model, u)
#{

#  # extract important initial quantities
#  beta <- model$coef
#  U <- beta %o% beta
#  M <- model$fisher
#  avar <- solve(M, symmetric = TRUE)
#  root <- model$root
#  pred <- model$pred
#  fam <- model$fam
#  x <- model$x
#  n <- nrow(x)
#  nnode <- ncol(x)
#  r <- length(beta)
#  modelmat <- matrix(model$modmat, ncol = r)

#  # obtain Gamma and efficiency gains
#  foo <- oneD(M = avar, U = U, u = u)
#  Gamma <- foo$G
#  mat <- foo$mat
#  ratio <- sqrt(diag(avar)/diag(mat))

#  # obtain the envelope estimator via profile likelihood
#  modmat.mnew <- modelmat %*% Gamma
#  arra <- array(modmat.mnew, c(n,nnode,u))
#  m2 <- aster(x, root, pred, fam, arra, type = "unconditional")
#  beta.env <- Gamma %*% m2$coef 

#  # see how many components of the envelope estimator reside
#  # within two SEs of the MLE
#  twoSE <- 2*sqrt(diag(avar)/n)
#  lower <- beta - twoSE  
#  upper <- beta + twoSE
#  cond <- (beta.env > lower) & (beta.env < upper)
#  cond.sat <- ifelse(cond,1,0)

#  # calculate a LRT  
#  env <-mlogl(beta.env, pred, fam, x, root, modmat = model$modmat,
#    type = "unconditional")$value
#  full <- mlogl(beta, pred, fam, x, root, modmat = model$modmat,
#    type = "unconditional")$value
#  LRT <- 2*(env - full)
#  pval <- pchisq(LRT, df = r - u, lower = F)

#  # calculate both AIC and BIC
#  k <- u*(r-u) + u + u*(u+1)/2 + (r-u)*(r-u+1)/2
#  k.full <- r + r*(r+1)/2
#  bic.env <- 2*env + k*log(n) 
#  aic.env <- 2*k + 2*env
#  bic.full <- 2*full + k.full*log(n)
#  aic.full <- 2*k.full + 2*full  
#  aic <- bic <- "full"
#  if(bic.env < bic.full) bic <- "envelope"
#  if(aic.env < aic.full) aic <- "envelope"
#  values <- matrix(c(aic.env, bic.env, aic.full, bic.full))
#  rownames(values) <- c("aic.env","bic.env","aic.full",
#    "bic.full")

#  # output
#  table <- cbind(beta, beta.env, cond.sat , ratio)
#  colnames(table) <- c("MLE", "envlp est.", "w/n 2 SE", 
#    "ratio")
#  out = list(eta = m2$coef, table = table, pval = pval, 
#    aic = aic, bic = bic, values = values)
#  return(out) 
#}
##########################################################




