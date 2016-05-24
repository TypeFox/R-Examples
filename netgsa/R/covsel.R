covsel <-
function(X, zero=NULL, one=NULL, lambda, rho=0.01, verbose=FALSE, eps = 1e-08) {
  n = dim(X)[1]
  p = dim(X)[2]
  Adj = matrix(0, p, p )
  
  if(is.null(zero)){
  	zero = matrix(0, p, p)
  }
  
  if (is.null(one)){
  	one = matrix(0, p, p)
  }
  
  ## To get the empirical covariance matrix
  X = scale(X, center = TRUE, scale = TRUE)
  for (i in 1:p) {
    Y = matrix(X[, i], ncol=1)
    Xmat = X[, -i]
    
    ## Get the zero and one indices.
    infoInd = one[i, -i] - zero[i, -i]
    beta = matrix(0, p - 1, length(lambda))
    
    if (sum((infoInd == 0)) == 0) {
      if (verbose) {
      	cat("Complete information known! \n")
      }
      beta[which(infoInd == 1), ] = 1
      beta[which(infoInd == -1), ] = 0
    } else {
      if (verbose) {
      	cat("Incomplete information known! \n")
      }
      
      if (sum((infoInd == -1)) == 0) { 
        if (sum((infoInd == 1)) == 0 ) {
          if (verbose) {
          	cat("Incomplete information: no 0's and no 1's! \n")
          }
          fit = glmnet(Xmat, Y, family = "gaussian", alpha = 1, lambda = lambda)
          beta = as.matrix(fit$beta)
        } else {
          if (verbose) {
          	cat("Incomplete information: no 0's, but with 1's! \n")
          }
          Xmat2 = Xmat[,(infoInd == 1)]
          Xmat1 = Xmat[,(infoInd != 1)]
          Pmat = Xmat2%*%pseudoinverse(t(Xmat2)%*%Xmat2)%*%t(Xmat2)
          Ynew = Y - Pmat%*%Y
          Xnew = Xmat1 - Pmat%*%Xmat1

          theta1 = as.matrix(glmnet(Xnew, Ynew, family = "gaussian", alpha = 1, lambda = lambda)$beta)
          theta2 = apply(theta1, 2,function(z){pseudoinverse(t(Xmat2)%*%Xmat2)%*%t(Xmat2)%*%(Y - Xmat1%*%z) })
          beta[(infoInd == 1),] = theta2
          beta[(infoInd != 1),] = theta1
        }
      } else {
        if (sum((infoInd == 1)) == 0 ) {
          if (verbose) {
          	cat("Incomplete information: with 0's and no 1's! \n")
          }
          beta[(infoInd == -1), ] = 0
          Xnew = Xmat[,(infoInd != -1)]

          fit = glmnet(Xnew, Y, family = "gaussian", alpha = 1, lambda = lambda)
          beta[(infoInd != -1),] = as.matrix(fit$beta)
        } else {
          if (verbose) {
          	cat("Incomplete information: with both 0's and 1's! \n")
          	} 
          beta[(infoInd==-1),] = 0
          Xmat2 = Xmat[,(infoInd == 1)]
          Xmat1 = Xmat[,(infoInd == 0)]
          Pmat = Xmat2%*%pseudoinverse(t(Xmat2)%*%Xmat2)%*%t(Xmat2)
          Ynew = Y - Pmat%*%Y
          Xnew = Xmat1 - Pmat%*%Xmat1

          theta1 = as.matrix(glmnet(Xnew, Ynew, family = "gaussian", alpha = 1, lambda = lambda)$beta)
          theta2 = apply(theta1, 2, function(z){pseudoinverse(t(Xmat2)%*%Xmat2)%*%t(Xmat2)%*%(Y - Xmat1%*%z) })
          beta[(infoInd == 1),] = theta2
          beta[(infoInd == 0),] = theta1
        }
      }
    }
    Adj[i, -i] = as.vector(beta)
  }
  
  ## Symmetrization
  Adj = (Adj + t(Adj))/2
  Adj = (abs(Adj) > eps)
  diag(Adj) = 0
  
  ## Estimate the weighted adjacency matrix based on graphical lasso
  info = zeroInd(Adj, 1)$zeroArr
  obj <- glasso(var(X), rho = rho, zero = info)
  
  sigma <- cov2cor(obj$w)
  siginv = chol2inv(chol(sigma))
  siginv[abs(siginv) < 1e-08] <- 0
  diag(siginv) <- 0
  
  return(list(Adj = Adj, wAdj = siginv))
}
