cv.covsel <-
function(X, zero=NULL, one=NULL, lambda, nfolds=NULL, verbose=FALSE) {
  n = dim(X)[1]
  p = dim(X)[2]
  Adj = array(0, c(p, p, length(lambda)))

  if(is.null(zero)){
  	zero = matrix(0, p, p)
  }
  
  if (is.null(one)){
  	one = matrix(0, p, p)
  }
  
  if (is.null(nfolds)){
  	nfolds = 5
  } 
    
  ## To get the empirical covariance matrix
  X = scale(X, center = TRUE, scale = TRUE)
  CVscore <- matrix(0, nfolds, length(lambda))
  cvIndex = cvFolds(n, K=nfolds,type="random")$which
  
  for (k in 1:length(unique(cvIndex))){
  	if (verbose){
  		cat("Fold",k, "\n")
  	}
  	trainD <- X[which(cvIndex!=k), ]
  	testD <- matrix(X[-which(cvIndex!=k), ], ncol=p) #could be a vector when it's leave-one-out-cv
  	
  	for (i in 1:p){
  		Y = matrix(trainD[, i], ncol=1)
  		Xmat = trainD[, -i]
  		infoInd = one[i, -i] - zero[i, -i]
	    beta = matrix(0, p - 1, length(lambda))
	    if (sum((infoInd == 0)) == 0) {
	        if (verbose ) {
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
            if (verbose ) {
            	cat("Incomplete information: no 0's, but with 1's! \n")
            }
            Xmat2 = Xmat[,(infoInd == 1)]
            Xmat1 = Xmat[,(infoInd != 1)]
            Pmat = Xmat2%*%pseudoinverse(t(Xmat2)%*%Xmat2)%*%t(Xmat2)
            Ynew = Y - Pmat%*%Y
            Xnew = Xmat1 - Pmat%*%Xmat1
            fit = glmnet(Xnew, Ynew, family = "gaussian", alpha = 1, lambda = lambda)
            
            theta1 = as.matrix(fit$beta)
            theta2 = apply(theta1, 2,function(z){pseudoinverse(t(Xmat2)%*%Xmat2)%*%t(Xmat2)%*%(Y - Xmat1%*%z) })
            beta[(infoInd == 1),] = theta2
            beta[(infoInd != 1),] = theta1
          }
        } else {
          if (sum((infoInd == 1)) == 0 ) {
            if (verbose ) {
            	cat("Incomplete information: with 0's and no 1's! \n")
            }
            beta[(infoInd == -1), ] = 0
            Xnew = Xmat[,(infoInd != -1)]
            fit = glmnet(Xnew, Y, family = "gaussian", alpha = 1, lambda = lambda)
            beta[(infoInd != -1),] = as.matrix(fit$beta)
          } else {
            if (verbose ) {
            	cat("Incomplete information: with both 0's and 1's! \n")
            } 
            beta[(infoInd==-1),] = 0
            Xmat2 = Xmat[,(infoInd == 1)]
            Xmat1 = Xmat[,(infoInd == 0)]
            Pmat = Xmat2%*%pseudoinverse(t(Xmat2)%*%Xmat2)%*%t(Xmat2)
            Ynew = Y - Pmat%*%Y
            Xnew = Xmat1 - Pmat%*%Xmat1
            fit = glmnet(Xnew, Ynew, family = "gaussian", alpha = 1, lambda = lambda)
            
            theta1 = as.matrix(fit$beta)
            theta2 = apply(theta1, 2, function(z){pseudoinverse(t(Xmat2)%*%Xmat2)%*%t(Xmat2)%*%(Y - Xmat1%*%z) })
            beta[(infoInd == 1),] = theta2
            beta[(infoInd == 0),] = theta1
          }
        }
      }
      resid = testD[, i] - testD[, -i]%*% as.matrix(beta)#a matrix of residuals
	  score.k = apply(resid, 2, function(t) sum(t^2))
      CVscore[k, ] = CVscore[k, ] + score.k      
  	}
  }
  CVscore <- apply(CVscore, 2, mean)
   
  return(list(lambda = fit$lambda, cve = CVscore))
}
