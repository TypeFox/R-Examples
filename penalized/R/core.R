###################################                                                           
# These functions are not user level functions
# They require a very specific input format and they 
# rely on the functions calling them for input checking
# These functions are not exported in the NAMESPACE file
###################################


###################################
# The core lasso and elastic net algorithm
###################################
.lasso <- function(beta, lambda, lambda2=0, positive, X, fit, trace = FALSE, 
  epsilon, maxiter) {

  # It is a general function for fitting L1-penalized models
  # possibly with an additional L2-penalty
  # Input:
  #  beta: a vector of length m (say) : starting values
  #  lambda: a vector of length m
  #  fit must be function(beta)
  #   Should return a list with at least:
  #     W:          The weights matrix, or its diagonal
  #     loglik:     The unpenalized loglikelihood, numeric)
  #     residuals:  The residuals
  
  m <- length(beta)
  n <- nrow(X)
  # are we fitting an elastic net?
  enet <- any(lambda2 != 0)   
  # find regression coefficients free of L1-penalty or positivity restraint
  free <- lambda == 0 & !positive      
  
  # initialize
  LL <- -Inf
  penalty <- penalty1 <- penalty2 <- Inf

  active <- !logical(m)
  nvar <- m

  tryNR <- FALSE
  NRfailed <- FALSE   
  whereNR <- NULL

  finished <- FALSE
  newfit <- TRUE
  retain <- 0.05
  cumsteps <- 0
  iter <- 0

  # iterate
  if (trace) cat("# nonzero coefficients:", m)
  while (!finished) {
    nzb <- (beta != 0)
    # calculate the local likelihood fit
    if (newfit) {
      activeX <- X[,nzb, drop=FALSE]
      linpred <- drop(activeX %*% beta[nzb])
      localfit <- fit(linpred)
      # Check for divergence
      if (is.na(localfit$loglik)) {
        if (trace) {
          cat(rep("\b", trunc(log10(nvar))+1), sep ="")
          warning("Model does not converge: please increase lambda.", call.=FALSE)
        }
        converged <- FALSE
        break
      }
      grad <- drop(crossprod(X, localfit$residuals))
      if (enet) {
        grad[active] <- grad[active] - lambda2[active] * beta[active]
      }
      oldLL <- LL
      oldpenalty <- penalty
      LL <- localfit$loglik
      penalty1 <- sum(lambda[active] * abs(beta[active]))
      if (enet) {
        penalty2 <- sum(lambda2[active] * beta[active] * beta[active])
      } else {
        penalty2 <- 0
      }
      penalty <- penalty1 + penalty2
      finishedLL <- (2 * abs(LL - oldLL) / (2 * abs(LL - penalty) + 0.1) < epsilon)
      finishedpen <- (2 * abs(penalty - oldpenalty) / (2 * abs(LL - penalty) + 0.1) < epsilon)
      cumsteps <- 0
    }

    # Calculate the penalized gradient from the likelihood gradient
    direction <- numeric(m)
    direction[nzb] <- grad[nzb] - lambda[nzb] * sign(beta[nzb])
    newb <- (!nzb) & ifelse(positive, grad > lambda, abs(grad) > lambda) 
    direction[newb] <- grad[newb] - lambda[newb] * sign(grad[newb])
    oldactive <- active   
    active <- nzb | newb
    activebeta <- beta[active]
    activedir <- direction[active]

    # check if retaining the old fit of the model does more harm than good
    oldnvar <- nvar
    nvar <- sum(active)
    if ((oldLL - oldpenalty > LL - penalty) || (nvar > 1.1* oldnvar)) {
      retain <- 0.5 * retain
    }

    # check convergence
    finishednvar <- !any(xor(active, oldactive))
    finished <- (finishedLL && finishedpen && finishednvar) || (all(activedir == 0)) || (iter == maxiter)

    if (!finished) {
      iter <- iter+1
      
      # Try Newton-Raphson, using the ridge routine if an L2 penalty is present
      if (tryNR) {
        activeX <- X[,active,drop=FALSE]
        if (enet && nvar > (n+1+sum(free)) ) {
          if (is.null(whereNR) || any(xor(whereNR, active))) {
            whereNR <- active
            P <- .makeP(activeX, lambda2[active], lambda[active], sign(activebeta))
            gams <- .solve(crossprod(t(P)), P %*% activebeta)
            PX <- P %*% t(activeX)
            Pl <- P * matrix(sqrt(lambda2[active]), nrow(P), ncol(P), byrow = TRUE)
            PlP <- crossprod(t(Pl))
          }
          if (is.list(localfit$W)) {
            PXdW <- PX * matrix(sqrt(localfit$W$diagW), nrow(PX), ncol(PX), byrow=TRUE)
            hessian <- - crossprod(t(PXdW)) + crossprod(t(PX %*% localfit$W$P)) - PlP
          } else if (length(localfit$W) > 1) {
            PXW <- PX * matrix(sqrt(localfit$W), nrow(PX), ncol(PX), byrow=TRUE)
            hessian <- -crossprod(t(PXW)) - PlP
          } else {
            hessian <- -crossprod(t(PX)) - PlP
          }
          Pgrad <- P %*% direction[active]
          shg <- drop(.solve(hessian, Pgrad))
          gams <- gams - shg
          NRbeta <- drop(crossprod(P, gams))
        } else { 
          if (is.list(localfit$W)) {
            XdW <- activeX * matrix(sqrt(localfit$W$diagW), nrow(activeX), ncol(activeX))
            hessian <- -crossprod(XdW) + crossprod(crossprod(localfit$W$P, activeX))
          } else if (length(localfit$W) > 1) {
            XW <- activeX * matrix(sqrt(localfit$W), nrow(activeX), ncol(activeX))
            hessian <- -crossprod(XW)
          } else {
            hessian <- -crossprod(activeX)
          } 
          if (enet) diag(hessian) <- diag(hessian) - lambda2[active]
          NRbeta <- activebeta - drop(.solve(hessian, direction[active]))
        }   
        NRfailed <- !all(sign(NRbeta) == sign(activebeta))
        if (!NRfailed) { 
          beta[active] <- NRbeta
          newfit <- TRUE
        } 
      } 

      if (!tryNR || NRfailed) {
        # find the second derivative of the likelihood in the projected direction
        if (newfit) {
          Xdir <- drop(X[,active, drop=F] %*% activedir)
          if (is.list(localfit$W)) {
            curve <- (sum(Xdir * Xdir * localfit$W$diagW) - drop(crossprod(crossprod(localfit$W$P, Xdir)))) / sum(activedir * activedir)
          } else if (length(localfit$W) > 1) {
            curve <- sum(Xdir * Xdir * localfit$W) / sum(activedir * activedir)
          } else {
            curve <- sum(Xdir * Xdir) / sum(activedir * activedir)
          }
          if (enet) {
            curve <- curve + sum(lambda2[active] * activedir * activedir) / sum(activedir * activedir)
          }
          topt <- 1 / curve
        }

        # how far can we go in the calculated direction before finding a new zero?
        tedge <- numeric(m)
        tedge[active] <- -activebeta / activedir
        tedge[tedge <= 0] <- 2 * topt
        tedge[free] <- 2* topt
        mintedge <- min(tedge)

        # recalculate beta
        if (mintedge + cumsteps < topt) {
          beta[active] <- activebeta + mintedge * activedir   
          beta[tedge == mintedge] <- 0  # avoids round-off errors
          cumsteps <- cumsteps + mintedge
          newfit <- (cumsteps > retain * topt) || (nvar == 1)  
          NRfailed <- FALSE
          tryNR <- FALSE
        } else {
          beta[active] <- activebeta + (topt - cumsteps) * activedir
          tryNR <- (cumsteps == 0) && !NRfailed && finishednvar && (enet || nvar < n)
          newfit <- TRUE
        }
      }
    } else {
      converged <- (iter < maxiter)
    }
    if (trace) {
      cat(rep("\b", max(1,trunc(log10(oldnvar))+1)), sep ="")
      cat(nvar)
      flush.console()
    }  
  }
  
  return(list(beta = beta, fit = localfit, penalty = c(L1 = penalty1, L2 = penalty2), iterations = iter, converged = converged))
}



###################################
# Adjusted lasso algorithm
# Tries to prevent large models by first fitting at higher values of lambda         
# Often faster for "pure" lasso
# Not recommended for elastic net
###################################
.steplasso <- function(beta, lambda, lambda2=0, positive, X, fit, trace = FALSE, epsilon, maxiter = Inf) {

  n <- nrow(X)
  finished <- FALSE
  while (!finished) {
    nzb <- (beta != 0)
    lp <- X[,nzb,drop=FALSE] %*% beta[nzb]
    gradient <- drop(crossprod(X[,!nzb,drop=FALSE], fit(lp)$residuals))
    rel <- gradient / lambda[!nzb]
    rel <- rel[rel>0 | !positive[!nzb]]
    if (length(rel) > n) {
      nextlambda <- sort(abs(rel), decreasing = TRUE)[n]
    } else {
      nextlambda <- 1
    }
    if (nextlambda <= 1) {
      finished <- TRUE
    }
    if(!finished) 
      out <- .lasso(beta, nextlambda * lambda, lambda2, positive, X, fit, trace, 1e-4, maxiter)
    else
      out <- .lasso(beta, lambda, lambda2, positive, X, fit, trace, epsilon, maxiter)
    beta <- out$beta
    if (trace && ! finished) cat(rep("\b", 24 + max(1,trunc(log10(sum(beta!=0))+1))), sep = "")
  }
  out
}


###################################
# Park & Hastie's suggestion for the next lambda where the active set changes
###################################
.park <- function(beta, lambda, lambda2=0, positive, X, fit) {

  gradient <- drop(crossprod(X[,drop=FALSE], fit$residuals))
  active1 <- ifelse(positive, (gradient - lambda) / lambda > -1e-4, (abs(gradient) - lambda) / lambda > -1e-4)
  active2 <- ifelse(positive, beta > 0, beta != 0)
  active <- active2 | active1
  if (sum(active) > 0) {
    activeX <- X[,active,drop=FALSE]
    if (is.list(fit$W)) {
      XadW <- activeX * matrix(sqrt(fit$W$diagW), nrow(activeX), ncol(activeX))  
      XaWXa <- crossprod(XadW) - crossprod(crossprod(fit$W$P, activeX))
      XdW <- X * matrix(sqrt(fit$W$diagW), nrow(X), ncol(X))  
      XWXa <- crossprod(XdW, XadW) - crossprod(crossprod(fit$W$P, X), crossprod(fit$W$P, activeX))    
    } else if (length(fit$W) > 1) {
      XaW <- activeX * matrix(sqrt(fit$W), nrow(activeX), ncol(activeX))
      XaWXa <- crossprod(XaW)
      XW <- X * matrix(sqrt(fit$W), nrow(X), ncol(X))
      XWXa <- crossprod(XW, XaW)
    } else {
      XaWXa <- crossprod(activeX)
      XWXa <- crossprod(X, activeX)
    } 
    signbetaactive <- ifelse(sign(beta[active]) == 0, sign(gradient[active]), sign(beta[active]))
    temp <- solve(XaWXa, signbetaactive * lambda[active])
    aa <- XWXa %*% temp / lambda
    hh1 <- (1 - gradient / lambda) / (1-aa) 
    hh1[hh1 < 0] <- Inf  
    hh2 <- (1 + gradient / lambda) / (1+aa)
    hh2[hh2 < 0] <- Inf
    hh3 <- - beta[active] / temp
    hh <- numeric(length(beta))
    hh[!active] <- pmin(hh1[!active], hh2[!active])
    hh[active] <- hh3      
    if (all(hh <= 0)) hh <- Inf else hh <- min(hh[hh>=0])
    if (hh < 1e-4) {hh <- 1e-4 }
  } else {
    hh <- 1e-4
  }
  newbeta <- beta
  newbeta[active] <- newbeta[active] + hh * temp
  list(hh = hh, beta = newbeta)
}


  
  
###################################
# The core ridge algorithm
###################################
.ridge <- function(beta, eta, Lambda, X, fit, trace = FALSE, epsilon = 1e-8, maxiter = 25) {

  if (missing(eta)) eta <- drop(X %*% beta)
  localfit <- fit(eta)
  
  iter <- 0
  oldLL <- -Inf
  finished <- FALSE

  while (!finished)
  {
    iter <- iter + 1
    if (trace) {
      cat(iter)
      flush.console()
    }
    if (is.matrix(Lambda)) {
      grad <- crossprod(X, localfit$residuals) - Lambda %*% beta
    } else {
      grad <- crossprod(X, localfit$residuals) - Lambda * beta
    }
    if (is.list(localfit$W)) {
      XdW <- X * matrix(sqrt(localfit$W$diagW), nrow(X), ncol(X))
      Hess <-  -crossprod(XdW) + crossprod(crossprod(localfit$W$P, X))
    } else if (length(localfit$W) > 1) {
      XW <- X * matrix(sqrt(localfit$W), nrow(X), ncol(X))
      Hess <- -crossprod(XW) 
    } else {  
      Hess <- -crossprod(X)
    }
    if (is.matrix(Lambda)) {
      Hess <- Hess - Lambda
    } else {
      diag(Hess) <- diag(Hess) - Lambda
    }
    shg <- drop(.solve(Hess, grad))
    beta <- beta - shg       
    eta <- drop(X %*% beta)

    if (trace) cat(rep("\b", trunc(log10(iter))+1), sep ="")
    
    localfit <- fit(eta)
    if (is.na(localfit$loglik) || iter == maxiter) {
      if (trace) {
        if (iter > 0) cat(rep("\b", trunc(log10(iter))+1), sep ="")
        warning("Model does not converge: please increase lambda.", call.=FALSE)
      }
      break
    } 
    if (is.matrix(Lambda)) {
      penalty <- as.numeric(0.5 * sum(beta * (Lambda %*% beta)))
    } else {
      penalty <- as.numeric(0.5 * sum(Lambda * beta * beta))
    }
    LL <- localfit$loglik - penalty
    
    # Check convergence
    finished <- ( 2 * abs(LL - oldLL) / (2 * abs(LL) + 0.1) < epsilon ) 
    half.step <- (LL < oldLL)
    if (half.step) {
      beta <- beta + 0.5 * shg
      eta <- drop(X %*% beta)
      localfit <- fit(eta)
      if (is.matrix(Lambda)) {
        penalty <- as.numeric(0.5 * sum(beta * (Lambda %*% beta)))
      } else {
        penalty <- as.numeric(0.5 * sum(Lambda * beta * beta))
      }
      LL <- localfit$loglik - penalty
    }
    oldLL <- LL
  }
  

  return(list(beta = beta, penalty = c(L1 = 0, L2 = penalty), fit = localfit, iterations = iter, converged = finished))
}


###################################
# Workhorse function for cross-validated likelihood
###################################
.cvl <- function(X, lambda1, lambda2, positive, beta, fit, chr, fusedl, groups, trace = FALSE, 
  betas = NULL, quit.if.failed = TRUE, save.predictions = TRUE, ...)  {

  n <- nrow(X)
  m <- ncol(X)

  # find the right fitting procedure
  useP <- FALSE
  
  if(!fusedl){
   if (all(lambda1 == 0) && !any(positive)) {
    if (m <= n) {
      cvfit <- function(leftout, beta) {
        subfit <- function(lp) fit$fit(lp, leftout)
        .ridge(beta = beta, Lambda = lambda2, X = X[!leftout,,drop = FALSE], fit = subfit, ...)
      }
    } else {
      useP <- TRUE
      P <- .makeP(X, lambda2)
      PX <- P %*% t(X)
      Pl <- P * matrix(sqrt(lambda2), nrow(P), ncol(P), byrow = TRUE)
      cvfit <- function(leftout, beta) {
        subfit <- function(lp) fit$fit(lp, leftout)
        gams <- .solve(tcrossprod(P), P %*% beta)
        PlP <- tcrossprod(Pl)
        out <- .ridge(beta = gams, Lambda = PlP, X = t(PX[,!leftout,drop = FALSE]), 
          fit = subfit, ...)
        out$beta <- drop(crossprod(P, out$beta))
        names(out$beta) <- names(beta)
        out
      }
    } 
  } else if (all(lambda2 == 0)) {
    cvfit <- function(leftout, beta) {
      subfit <- function(lp) fit$fit(lp, leftout)
      .steplasso(beta = beta, lambda = lambda1, positive = positive, 
        X = X[!leftout,,drop = FALSE], fit = subfit, ...)
    }
  } else {
    cvfit <- function(leftout, beta) {
      subfit <- function(lp) fit$fit(lp, leftout)
      .lasso(beta = beta, lambda = lambda1, lambda2 = lambda2, positive = positive, 
        X = X[!leftout,,drop = FALSE], fit = subfit, ...)
    }
   }
  } else if (fusedl){
    cvfit <- function(leftout, beta) {
      subfit <- function(lp) fit$fit(lp, leftout)
      .flasso(beta = beta, lambda1 = lambda1, chr = chr, lambda2 = lambda2, positive = positive,
        X = X[!leftout,,drop = FALSE], fit = subfit,...)
    }
  }

  # "groups" input lists fold allocation for each subject %in% 1:fold
  fold <- max(groups)
                                      
  # fit the full model and make an m x fold matrix of beta if necessary
  fullfit <- cvfit(logical(n), beta)
  if (is.null(betas)) {
    if (fullfit$converged) 
      betas <- matrix(fullfit$beta, m, fold)
    else
      betas <- matrix(beta, m, fold)
  } 
  
  # True cross-validation starts here
  if (fold > 1) {
    failed <- FALSE
    predictions <- vector("list", fold)
    cvls <- sapply(1:fold, function(i) {
      if (!failed) {
        if (trace) {
          cat(i)
          flush.console()
        }
        leaveout <- (groups == i)
        foldfit <- cvfit(leaveout, betas[,i])
        lin.pred <- numeric(n)
        lin.pred[leaveout] <- X[leaveout, foldfit$beta != 0, drop=FALSE] %*% 
          foldfit$beta[foldfit$beta != 0]
        lin.pred[!leaveout] <- foldfit$fit$lp0
        if (save.predictions)
          predictions[[i]] <<- fit$prediction(lin.pred[leaveout], foldfit$fit$nuisance, leaveout)
        betas[,i] <<- foldfit$beta
        if (trace) cat(rep("\b", trunc(log10(i))+1), sep ="")
        out <- fit$cvl(lin.pred, leaveout)
        if (quit.if.failed && (is.na(out) || abs(out) == Inf || foldfit$converged == FALSE)) failed <<- TRUE
      } else {
        out <- NA
      }
      
      out
    })
                                      
    if (failed || any(is.na(cvls))) cvls <- -Inf
    
  } else {
    cvls <- NA
    predictions <- NA
    betas <- NA
  }

  list(cvl = sum(cvls), cvls=cvls, fit = fullfit, betas = betas, predictions = predictions)
}

#######################################
# makes a reduced basis for L2-penalized newton-raphson in the p > n case
#######################################
.makeP <- function(X, lambda2, lambda1 = 0, signbeta) {

  n <- nrow(X)
  p <- ncol(X)
                                                      
  free2 <- (lambda2 == 0)
  free1 <- all(lambda1 == 0)
  m <- sum(free2)
  
  if (free1) {
    P <- matrix(0, n+m, p)
  } else {
    P <- matrix(0, n+m+1, p)
  }

  # First columns: free variables in case of no L2-penalization
  for (i in seq_len(m)) P[i, which(free2)[i]] <- 1
                                                
  # Next n columns: column span of X
  P[m + 1:n, which(!free2)] <- X[,!free2,drop=FALSE] * matrix(1/lambda2[!free2], n, p-m, byrow=TRUE)

  # Additional column due to L1-penalization
  if (!free1) {          
    P[n+m+1,which(!free2)] <- (lambda1*signbeta/lambda2)[!free2]
  }
                    
  # check for singularity
  rownames(P) <- paste("V", 1:nrow(P))
  qrP <- qr(t(P))  # not efficient?
  if (qrP$rank < nrow(P)) {
    qR <- qr.R(qrP)
    keep <- colnames(qR)[sort.list(abs(diag(qR)), decreasing=TRUE)[1:qrP$rank]]
    P <- P[keep,]
  }
  colnames(P) <- NULL
  
  # Numerical stabilization          
  #P <- P / matrix(apply(P, 1, sd), nrow(P), ncol(P), byrow=F)

  return(P)
}

#######################################
# a solve() function that does not complain about near-singularity
# often dangerous, but very useful here
#######################################
.solve <- function(a,b) {
  out <- try(qr.coef(qr(a, LAPACK=TRUE), b))
  if (is(out, "try-error")) stop("Matrix inversion failed. Please increase lambda1 and/or lambda2", call. = FALSE)
  return(out)
}

#######################################
# calculate cross-validation folds
#######################################
.getFolds <- function(fold, n) {
  if(length(fold) == 1) {
    if (fold == n) {
      groups <- 1:n
    } else {
      groups <- sample(n) %% fold + 1
    }
  } else {
    if (length(fold) == n) {
      groups <- fold 
      fold <- max(groups)
    } else {
      stop("incorrect input of \"fold\"", call.=FALSE)
    }
  }
  if (!all(1:fold %in% groups)) stop("incorrect input of \"fold\"", call.=FALSE)
  return(groups)
}

.cvlapprox <- function(X, lambda1, lambda2, positive, beta, fit, groups,trace=FALSE,
  betas=NULL, quit.if.failed=TRUE, save.predictions=TRUE, ...)
{
  n <- nrow(X)
  m <- ncol(X)

  fold <- max(groups)

  if(fold==n)   #if loocv, make sure vector "groups" is ordered
  {
    groups <- c(1:n)
  }

  useP <- FALSE
  if(m<=n)  #use beta
  {
      fullfit = .ridge(beta = beta, Lambda = lambda2, X = X, fit = fit$fit)
            
      if(!fullfit$converged)
      {
        cvls <- -Inf
        predictions <- NA
        betas <- NA
        return(list(cvl = sum(cvls), cvls=cvls, fit = fullfit, betas = betas, predictions = predictions))  
      }
     
      beta <- fullfit$beta
      linpreds <- fullfit$fit$lp
      delta <- fullfit$fit$residuals
      W <- fullfit$fit$W
      
      #check corresponding model: logistic/poisson, linear or Cox by inspection of W (either diag, 1 or list)

      if (is.list(W)) #cox regression
      {
          diagW = W$diagW  
          #W = diag(diagW)
          beta_int <- c(0,beta)
          Xint <- cbind(1,X)
          # XWX = crossprod(W^{1/2}X)
          # W^{1/2}X = matrix(sqrt(dW), nrow(X), ncol(X)) * X
          # inv <- solve(t(Xint)%*%W%*%Xint+diag(c(0,lambda))) 
          # inv <- solve(crossprod(matrix(sqrt(diagW),nrow(X),ncol(X))*X)+diag(c(0,lambda)))
          XWXlambda <- crossprod(matrix(sqrt(diagW),nrow(Xint),ncol(Xint))*Xint)
          diag(XWXlambda) <- diag(XWXlambda) + c(0,lambda2)
          inv_tXint <- solve(XWXlambda,t(Xint))
          WXint <- matrix(diagW,nrow(Xint),ncol(Xint))*Xint
          
          if(fold < n) #kfold
          {
              betas <- matrix(0,(m+1),fold)
              
              for(i in 1:fold)
              {
                  #beta_k <- .getBeta(Xint,beta_int,groups,i,W,delta,inv)
                  beta_k <- .getBeta2(WXint,beta_int,groups,i,delta,inv_tXint)
                  betas[,i] <- beta_k
              }          
          }
          else #loocv
          {
              #V <- W^(1/2)%*%Xint%*%inv%*%t(Xint)%*%W^(1/2) 
              V <- WXint%*%inv_tXint  
              betas <- matrix(beta_int,length(beta_int),fold)
              #betas_int <- betas_int - inv%*%t(Xint)%*%diag(1/(1-diag(V))*delta)
              betas <- betas - inv_tXint*matrix((1/(1-diag(V))*delta),nrow(inv_tXint),ncol(inv_tXint),byrow=TRUE) 
          }
          X <- Xint   # nodig voor lin.pred
      } 
      else if (length(W) > 1)  #poisson/logistic regression 
      {
          diagW <- W                                                            #NB: W is just the diagonal
          XWXlambda <- crossprod(matrix(sqrt(diagW),nrow(X),ncol(X))*X)       
          diag(XWXlambda) <- diag(XWXlambda) + lambda2
          inv_tX <- solve(XWXlambda,t(X))
          WX <- matrix(diagW,nrow(X),ncol(X))*X                               
          #inv <- solve(t(X)%*%W%*%X+diag(lambda))

          if(fold < n) #kfold
          {
              betas <- matrix(0,m,fold)
              
              for(i in 1:fold)
              {
                  beta_k <- .getBeta2(WX,beta,groups,i,delta,inv_tX) 
                  betas[,i] <- beta_k
              }          
          }
          else #loocv
          {
              V <- WX%*%inv_tX  #forceSymmetric could be useful
              betas <- matrix(beta,m,fold)
              betas <- betas - inv_tX*matrix((1/(1-diag(V))*delta),nrow(inv_tX),ncol(inv_tX),byrow=TRUE) 
      
          }
      } 
      else  #linear regression 
      {
          XXlambda <- crossprod(X)
          diag(XXlambda) <- diag(XXlambda) + lambda2
          inv_tX <- solve(XXlambda,t(X))
          
          #inv <- solve(t(X)%*%X+diag(lambda))
          
          if(fold < n) #kfold
          {
              betas <- matrix(0,m,fold)
              
              for(i in 1:fold)
              {
                beta_k <- .getBetaLinear2(X,beta,groups,i,delta,inv_tX)
                betas[,i] <- beta_k
              }          
          }
          else #loocv
          {
              V <- X%*%inv_tX  
              betas <- matrix(beta,m,fold)
              betas <- betas - inv_tX*matrix((1/(1-diag(V))*delta),nrow(inv_tX),ncol(inv_tX),byrow=TRUE) 
          }
      }

  }
  else #use gamma
  {
      useP <- TRUE
      P <- .makeP(X, lambda2)
      PX <- P %*% t(X)
      Pl <- P * matrix(sqrt(lambda2), nrow(P), ncol(P), byrow = TRUE)
      gams <- .solve(tcrossprod(P), P %*% beta)  #initial gamma 
      PlP <- tcrossprod(Pl)
      
      fullfit <- .ridge(beta = gams, Lambda = PlP, X = t(PX), fit = fit$fit)
      gams <- fullfit$beta    #new gamma
      fullfit$beta <- drop(crossprod(P, gams))  #prepare fullfit for return statement
      names(fullfit$beta) <- names(beta)
      
      if(!fullfit$converged)
      {
        cvls <- -Inf
        predictions <- NA
        betas <- NA
        return(list(cvl = sum(cvls), cvls=cvls, fit = fullfit, betas = betas, predictions = predictions)) 
      }      
      
      linpreds <- fullfit$fit$lp   
      delta <- fullfit$fit$residuals
      W <- fullfit$fit$W
      
      #check which model you're in

      if (is.list(W)) #cox regression
      {
          diagW <- W$diagW  #vector
          G <- P
          k <- nrow(G)
          l <- ncol(G)

          G <- cbind(c(1,rep(0,k)),rbind(rep(0,l),G))
          tG <- t(G)
          gamma_int <- c(0,gams)
          Xint <- cbind(m,X) #NB, m ipv 1!      m=median(X)
          #B <- Xint%*%tG
          B <- tcrossprod(Xint,G)
          #Atilde <- G%*%diag(c(0,lambda2))%*%tG 
          Gl <- G * matrix(sqrt(c(0,lambda2)), nrow(G), ncol(G), byrow = TRUE)  
          Atilde <- tcrossprod(Gl)
          
          BWB <- crossprod(matrix(sqrt(diagW),nrow(B),ncol(B))*B)
          BWBlambda <- BWB + Atilde
          inv_tB <- solve(BWBlambda,t(B))
          WB <- matrix(diagW,nrow(B),ncol(B))*B
          
          #inv <- solve(t(B)%*%W%*%B+Atilde)
          
          if(fold < n) #kfold
          {
              betas <- matrix(0,(m+1),fold)   
              
              for(i in 1:fold)
              {
                gamma_k <- .getBeta2(WB,gamma_int,groups,i,delta,inv_tB)
                beta_k <- tG %*% gamma_k
                betas[,i] <- beta_k
              } 
          }
          else #loocv
          {
              V <- WB %*% inv_tB  #forceSymmetric?
              gammas_int <- matrix(gamma_int,length(gamma_int),fold)
              gammas_int <- gammas_int - inv_tB*matrix((1/(1-diag(V))*delta),nrow(inv_tB),ncol(inv_tB),byrow=TRUE) 
              betas <- tG %*% gammas_int
          }
          X <- Xint
      } 
      else if (length(W) > 1)  #poisson/logistic regression 
      {
          diagW <- W
          G <- P
          tG <- t(G)
          B <- t(PX)
          Atilde <- PlP
          BWB <- crossprod(matrix(sqrt(diagW),nrow(B),ncol(B))*B)
          BWBlambda <- BWB + Atilde
          inv_tB <- solve(BWBlambda,t(B))
          WB <- matrix(diagW,nrow(B),ncol(B))*B
          #inv <- solve(t(B)%*%W%*%B+Atilde)
          
          if(fold < n) #kfold
          {
              betas <- matrix(0,m,fold)
              
              for(i in 1:fold)
              {
                gamma_k <- .getBeta2(WB,gams,groups,i,delta,inv_tB)
                betas[,i] <- tG %*% gamma_k
              }          
          }
          else #loocv
          {
              V <- WB %*% inv_tB  #forceSymmetric?
              gammas <- matrix(gams,length(gams),fold)
              gammas <- gammas - inv_tB*matrix((1/(1-diag(V))*delta),nrow(inv_tB),ncol(inv_tB),byrow=TRUE) 
              betas <- tG %*% gammas      
          }
      } 
      else  #linear regression    
      {  
          G <- P
          tG <- t(G)
          B <- t(PX)
          Atilde <- PlP
         
          BBlambda <- crossprod(B) + Atilde
          inv_tB <- solve(BBlambda,t(B))
          #inv <- solve(PX%*%t(PX)+PlP)  
    
          if(fold < n) #kfold
          {
              betas <- matrix(0,m,fold)
              
              for(i in 1:fold)
              {
                  gamma_k <- .getBetaLinear2(B,gams,groups,i,delta,inv_tB)
                  betas[,i] <- tG %*% gamma_k
              }          
          }
          else #loocv
          {
              V <- B%*%inv_tB
              gammas <- matrix(gams,length(gams),fold)
              gammas <- gammas - inv_tB*matrix((1/(1-diag(V))*delta),nrow(inv_tB),ncol(inv_tB),byrow=TRUE)
              betas <- tG %*% gammas
          }
      }      
    
  
  }

  # True cross-validation starts here
  if (fold > 1) 
  {
    failed <- FALSE
    predictions <- vector("list", fold)
    cvls <- sapply(1:fold, function(i) 
    {
      if (!failed) 
      {
          leaveout <- (groups == i)
          lin.pred <- X %*% betas[,i] #? 

          if (save.predictions)
          {
            if (length(W)==1)   #linear model
            {
              y <- fullfit$fit$residuals + fullfit$fit$lp
              rss <- drop(crossprod(y[!leaveout]-X[!leaveout,]%*%betas[,i]))
              nr <- n - sum(leaveout)
              predictions[[i]] <<- fit$prediction(lin.pred[leaveout], list(sigma2 = rss/nr), leaveout)           
            } 
            else
            {
              predictions[[i]] <<- fit$prediction(lin.pred[leaveout], fullfit$fit$nuisance, leaveout) #in cox: lin.pred contains intercept            
            }
          }
          
          out <- fit$cvl(lin.pred, leaveout)
          
          if (quit.if.failed && (is.na(out) || abs(out) == Inf || fullfit$converged == FALSE)) failed <<- TRUE   
      } 
      else 
      {
          out <- NA
      }
      
      out
      
    })
                                      
    if (failed || any(is.na(cvls))) cvls <- -Inf
    
  } 
  else 
  {
    cvls <- NA
    predictions <- NA
    betas <- NA
  }

  list(cvl = sum(cvls), cvls=cvls, fit = fullfit, betas = betas, predictions = predictions)
}

#only difference with getBeta2 is the weightmatrix W that does not occur in the linear model
.getBetaLinear2 <- function(X,beta,groups,j,delta,inv_tX)
{
  Xj <- X[groups==j,,drop=FALSE]
  inv_tXj <- inv_tX[,groups==j,drop=FALSE]
  Imin <- -Xj%*%inv_tXj       #forceSymmetric
  diag(Imin) <- diag(Imin+1)

  betaj <- beta - inv_tXj%*%solve(Imin,delta[groups==j])      #Imin can be singular in some specific cases.. 
  #betaj <- beta - inv%*%tXj%*%solve(I-Xj%*%inv%*%tXj)%*%delta[groups==j]
  return(betaj)
}

.getBeta2 <- function(WX,beta,groups,j,delta,inv_tX) 
{
  WXj <- WX[groups==j,,drop=FALSE]
  inv_tXj <- inv_tX[,groups==j,,drop=FALSE]
  Imin <- -WXj%*%inv_tXj    #forceSymmetric
  diag(Imin) <- diag(Imin+1)

  betaj <- beta - inv_tXj%*%solve(Imin,delta[groups==j])
  return(betaj)
}

###Fused Lasso
####################

###############################################################################################
.flasso <- function(beta, chr, lambda1, lambda2,fit, X, positive, trace = FALSE,
  epsilon=1e-5, maxiter=Inf) {

  # Input:
  #  beta: a vector of length m (say) : starting values
  #  lambda1 and lambda2: vectors of length m
  #  X : data matrix
  #  positive: logical vector denoting the coefficients which should pe penalized
  #   Should return a list with at least:
  #     W:          The weights matrix, or its diagonal
  #     loglik:     The unpenalized loglikelihood, numeric)
  #     residuals:  The residuals
  # If the model includes an Intercept then starting with all the beta coefficients
  # equal to zero might slow down the convergence. Hence it's better to use a
  # warm start

  m=length(beta)
  n=nrow(X)
  

# find regression coefficients free of L1-penalty or positivity restraint
  
  free <- (lambda1 == 0 & lambda2==0) & !positive

# initialize
  
  LL <- -Inf
  penalty <- penalty1 <- penalty2 <- Inf
  converged = FALSE
  active <- !logical(m)

  if(any(free)){direc=numeric(length(free[-which(free)]))}
  if(!any(free)){direc=numeric(length(beta))}
  if(all(free)){direc=numeric(length(free))}

  actd = (diff(direc))!=0
  checkd = matrix(0,1,2)
  
  nvar <- m
  oldmo=length(beta)
  
  tryNR <- FALSE
  NRfailed <- FALSE
  finish <- FALSE
  newfit <- TRUE
  
  retain <-0
  cumsteps <- 0
  iter <- 0
 #iterate

while(!finish){
        
           nzb = (beta!=0)
            if(any(free)){direc=numeric(length(free[-which(free)]))}
            if(!any(free)){direc=numeric(length(beta))}
            if(all(free)){direc=numeric(length(free))}
    # calculate the local likelihood fit
  
if (newfit) {
        activeX <- X[,nzb, drop=FALSE]
        linpred <- drop(activeX %*% beta[nzb])
        lp=linpred
        localfit=fit(lp)
 # Check for divergence
      if (is.na(localfit$loglik)) {
        if (trace) {
          cat(rep("\b", trunc(log10(nvar))+1), sep ="")
          warning("Model does not converge: please increase lambda.", call.=FALSE)
        }
        converged <- FALSE
        break
      }
    grad <- drop(crossprod(X,localfit$residuals))
    oldLL <- LL
    oldpenalty <- penalty
    LL <- localfit$loglik
    penalty1=sum(lambda1[active]*abs(beta[active]))
    penalty2 = 0
  
if((!any(chr==0) && length(table(chr))<=1) || (any(chr ==0) && length(table(chr))<=2)){
      penalty2= sum(lambda2[1:(length(lambda2)-1)]*abs(((diff(beta)))))
      penalty=penalty1+penalty2
      finishedLL <- (2 * abs(LL - oldLL) / (2 * abs(LL - penalty) + 0.1) < epsilon)
      finishedpen <- (2 * abs(penalty - oldpenalty) / (2 * abs(LL - penalty) + 0.1) < epsilon)
      cumsteps <- 0
    }

if((!any(chr==0) && length(table(chr))>1) || (any(chr==0) && length(table(chr))>2)){
     chrn=chr[!free]
      
     for (chri in 1:length((as.numeric(names(table(chrn)))))){
      betac=beta[!free]
      lambda2c=lambda2[!free]
      beta1.1 = betac[which(chrn==(as.numeric(names(table(chrn))))[chri])]
      lambda2.1 = lambda2c[which(chrn==(as.numeric(names(table(chrn))))[chri])]
      penalty21= sum(lambda2.1[1:(length(lambda2.1)-1)]*abs(((diff(beta1.1)))))
      penalty2= penalty2 + penalty21
      }
  
      penalty <- penalty1 + penalty2
      finishedLL <- (2 * abs(LL - oldLL) / (2 * abs(LL - penalty) + 0.1) < epsilon)
      finishedpen <- (2 * abs(penalty - oldpenalty) / (2 * abs(LL - penalty) + 0.1) < epsilon)
      cumsteps <- 0
    }
}



    
# Calculate the penalized gradient from the likelihood gradient
    ###Initialize###
   
     

    #####  Core ########################
chrn=chr[!free]

if(!(all(free))){

       gradc=grad[!free]
       nzbc=nzb[!free]
       betac=beta[!free]
       lambda1c=lambda1[!free]
       lambda2c=lambda2[!free]
       positivec=positive[!free]

   for (chri in 1:length((as.numeric(names(table(chrn)))))){

     beta.1 = betac[which(chrn==(as.numeric(names(table(chrn))))[chri])]
     grad.1 = gradc[which(chrn==(as.numeric(names(table(chrn))))[chri])]
     lambda1.1 = lambda1c[which(chrn==(as.numeric(names(table(chrn))))[chri])]
     lambda2.1 = lambda2c[which(chrn==(as.numeric(names(table(chrn))))[chri])]
     nzb.1 = nzbc[which(chrn==(as.numeric(names(table(chrn))))[chri])]
     positive.1 = positivec[which(chrn==(as.numeric(names(table(chrn))))[chri])]

     gdir <- getdirec(grad=grad.1,nzb=nzb.1,beta=beta.1,lambda1=lambda1.1,lambda2=lambda2.1,positive=positive.1)
     direc[which(chrn==(as.numeric(names(table(chrn))))[chri])] = gdir$direc_i

        }


     if(any(free)){direc=c(grad[which(free)],direc)}
}  else {direc=grad}



     oldcheckd = checkd
     oldactive=active
     active= direc!=0 | nzb

   
     checkd = nonZero(direc[!free])
     actd = (diff(direc[!free]))!= 0


     activdir=direc[active]
     activbeta=beta[active]


# check if retaining the old fit of the model does more harm than good
    oldnvar <- nvar
    nvar <- sum(active)
    if ((oldLL - oldpenalty > LL - penalty) || (nvar > 1.1* oldnvar)) {
      retain <- 0.5 * retain
    }
  # check convergence
    if(length(checkd)==length(oldcheckd)){
    finishednvar <- (!any(xor(active, oldactive)) &&  all(checkd==oldcheckd))
    }
    if(length(checkd)!=length(oldcheckd)){
    finishednvar <- (!any(xor(active, oldactive)) &&  length(checkd)==length(oldcheckd))
    }

    finishedn <- !any(xor(active, oldactive))
    if(any(free)){
    finish <- (finishedLL && finishedn && finishedpen) || (all(activdir == 0)) || (iter == maxiter)}
    if(!any(free)){
    finish <- (finishedLL && finishedn && finishedpen) || (all(activdir == 0)) || (iter == maxiter)}


  if (!finish) {
      iter <- iter+1
  
      
              if (tryNR){
             NRs<- tryNRstep(beta=beta,direc=direc,active=active,X=X,newfit=newfit,localfit=localfit)
             beta = NRs$beta
             newfit = NRs$newfit
             NRfailed = NRs$NRfailed
           }

        
               


 if (!tryNR || NRfailed) {
 
      if (newfit) {
          Xdir <- drop(X[,active,drop=F] %*% activdir)
          if (is.list(localfit$W)) {
            curve <- (sum(Xdir * Xdir * localfit$W$diagW) - drop(crossprod(crossprod(localfit$W$P, Xdir)))) / sum(activdir * activdir)
          } else if (length(localfit$W) > 1) {
            curve <- sum(Xdir * Xdir * localfit$W) / sum(activdir * activdir)
          } else {
            curve <- sum(Xdir * Xdir) / sum(activdir * activdir)
          }
           topt <- 1 / curve
       }
       
 

        # how far can we go in the calculated direction before finding a new zero?
        tedge <- numeric(length(beta))
        tedge[active] <- -activbeta / activdir
        tedge[tedge <= 0] <-tedge[free] <-  2 * topt
        if(length(lambda1)==1 && length(lambda2)==1 && lambda1==0 && lambda2==0){tedgeg=NULL}
        if((length(lambda1)!=1 || length(lambda2)!=1) && (sum(lambda1==0)!=length(lambda1) || sum(lambda2==0)!=length(lambda2))){
          tedgeg = -diff(beta)/diff(direc)
          tedgeg[which(is.na(tedgeg))]=0
          dchr = diff(chr)
          tedgeg[which(dchr!=0)]=0
          tedgeg[tedgeg <= 0] <-  2 * topt
          }


         if(!is.null(tedgeg)){
          mintedge <- min(min(tedge),min(tedgeg))}
          if(is.null(tedgeg)){mintedge <- min(tedge)}
    
      
        # recalculate beta
          if (mintedge + cumsteps < topt) {
          
          beta[active] <- activbeta + mintedge * activdir
          beta[tedge == mintedge] <- 0  # avoids round-off error
          if(!is.null(tedgeg)){
          tmin=which(tedgeg==mintedge)
          beta[tmin]=beta[tmin+1]
           }
          cumsteps=cumsteps+mintedge
          newfit <- (cumsteps > retain * topt) || (nvar == 1)  
          NRfailed <- FALSE
          tryNR <- FALSE

          } else {
          beta[active] <- activbeta + (topt - cumsteps) * activdir
          #if(iter<=1){tryNR <- (cumsteps == 0) && finishednvar && !NRfailed && nvar < m}else{tryNR <- (cumsteps == 0)&& !NRfailed && nvar < n}
          tryNR <- (cumsteps == 0) && !NRfailed && finishednvar && nvar < m
          newfit <- TRUE
          }
          
        }

      }
      
  else {
      converged <- (iter < maxiter)
      }

 
     
     }

    
  return(list(beta = beta, fit = localfit, penalty = c(L1 = penalty1, L2 = penalty2), iterations = iter, converged = converged))
}
######################################################################################
##function for getting the indices of non-zero coefficients and coefficients equal to each other#####
nonZero <- function(inp)
{
 n <- length(inp)

 indices <- matrix(rep(0,n*2),n,2)
 j=1
 i=1
 while(i<=n)
 {
  if(inp[i] == 0)
  {
    i <- i+1
  }
  else if(i<n && (inp[i+1]==inp[i]))
  {
    indices[j,1] <- i
    while(i<n && inp[i+1]==inp[i])
    {
      i <- i+1
    }
    indices[j,2] <- i
    j <- j+1
    i <- i+1
  }
  else
  {
    indices[j,1] <- indices[j,2] <- i
    i <- i+1
    j <- j+1
  }

 }

 return(indices[(1:j-1),])

}
#########################################


###################################

getdirec <- function(grad,nzb,beta,lambda1,lambda2,positive){



    #####################
    #####  Core ########################
    direc_i=numeric(length(beta))
     i=1
     vmax = 0
     lam1 <-  lam2 <- P_grad<-cu_g <-numeric(length(beta))
     oldvmax = vmax
     olddirec_i=direc_i
     if(any(nzb)){
    nz <- nonZero(beta)
              nz=rbind(nz,c(0,0))
              for (ik in 1:(nrow(nz)-1)){
              cu_g[nz[ik,1]:nz[ik,2]]=(sum(grad[nz[ik,1]:nz[ik,2]]))
              P_grad[nz[ik,1]:nz[ik,2]]=(cu_g[nz[ik,1]:nz[ik,2]])/(length(nz[ik,1]:nz[ik,2]))
               if(nz[ik,1]==nz[ik,2]){
                lam1[nz[ik,1]]=lambda1[nz[ik,1]]*sign(beta[nz[ik,1]])
               if(nz[ik,1]==1){
                   lam2[nz[ik,1]]=-lambda2[nz[ik,1]]*sign(beta[nz[ik,1]+1]-beta[nz[ik,1]])*-1
               }
               if(nz[ik,1]==length(beta)){
                   lam2[nz[ik,1]]=-lambda2[nz[ik,1]]*sign(beta[nz[ik,1]]-beta[nz[ik,1]-1])
               }
               if(nz[ik,1]!=1 && nz[ik,1]!=length(beta)){
               lam2[nz[ik,1]]=((-lambda2[nz[ik,1]]*sign(beta[nz[ik,1]]-beta[nz[ik,1]-1]))-(lambda2[nz[ik,1]]*sign(beta[nz[ik,1]+1]-beta[nz[ik,1]])*-1))
                              }
                 }
                 if(nz[ik,1]!=nz[ik,2]){
                  lam1[nz[ik,1]:nz[ik,2]]=lambda1[nz[ik,1]:nz[ik,2]]*sign(beta[nz[ik,1]:nz[ik,2]])
                   if(nz[ik,1]==1 && nz[ik,2]!=length(beta)){
                     lam2[nz[ik,1]:nz[ik,2]]=(-lambda2[nz[ik,1]]*sign(beta[nz[ik,2]+1]-beta[nz[ik,1]])*-1)/(length(nz[ik,1]:nz[ik,2]))
                     }
                   if(nz[ik,1]!=1 && nz[ik,2]==length(beta) ){
                   lam2[nz[ik,1]:nz[ik,2]]=(-lambda2[nz[ik,1]]*sign(beta[nz[ik,1]]-beta[nz[ik,1]-1]))/(length(nz[ik,1]:nz[ik,2]))
                   }
                   if(nz[ik,2]==length(beta) && nz[ik,1]==1 ){
                   lam2[nz[ik,1]:nz[ik,2]]=0
                   }
                   if(nz[ik,2]!=length(beta) & nz[ik,1]!=1){
                   lam2[nz[ik,1]:nz[ik,2]]=(((-lambda2[nz[ik,1]]*sign(beta[nz[ik,1]]-beta[nz[ik,1]-1]))-(lambda2[nz[ik,1]]*sign(beta[nz[ik,2]+1]-beta[nz[ik,1]])*-1)))/(length(nz[ik,1]:nz[ik,2]))
                    }
                  }
                 }
                 }
 direc_i= P_grad - lam1 + lam2
 olddirec_i=direc_i


  while(i<length(beta)){
      lam1_1 <- lam1_2  <- lam2_1 <- lam2_2 <-cu_g1 <- cu_g2 <- P_grad1 <- P_grad2<-numeric(length(beta))
       direc_i <- numeric(length(beta))
       direc_i= P_grad - lam1 + lam2

       

  ######## Penalized Gradient for nonzero coefficients##################
  if(i<length(beta)){

      if(nzb[i]){
       nzb_m=c(nzb[i:length(beta)],F)
       m=which(diff(nzb_m)==-1)
       m=i+(m-1)
       m=m[1]
                                  #############gives the length of the non-zero stretch######
      if(i==length(nzb)){m=length(nzb)}
      betam=beta[i:m]
        if(sum(diff(betam))==0){j=m}
        if(sum(diff(betam))!=0){
           j=which(diff(betam)!=0)
           j=i+(j-1)
           j=j[1]}
           m=j

    }

    if(!nzb[i]){

            nzb_mo=c(nzb[i:length(beta)],T)
            mo1=which(diff(nzb_mo)==1)[1]
            mo1=i+(mo1-1)
            m=mo1

            if(i==length(nzb)){m=length(beta)}

          }

      P_grad1=P_grad
      lam1_1=lam1
      lam2_1=lam2
      P_grad1[i:m]=0
      lam1_1[i:m]=0
      lam2_1[i:m]=0

if((length(i:m)==1) && (!nzb[i]) ){

    P_grad1[m]=grad[m]
       lam1_1[m]=lambda1[m]*sign(grad[m])
      if(m!=length(beta)){
        if(m!=1){
          if(beta[m-1]!=0 & beta[m+1]!=0){
            lam2_1[m]=((-lambda2[m]*sign(beta[m]-beta[m-1]))-(lambda2[m]*sign(beta[m+1]-beta[m])*-1))

           }
           if(beta[m-1]==0 & beta[m+1]!=0){
                lam2_1[m]=((-lambda2[m]*sign(grad[m]-0))-(lambda2[m]*sign(beta[m+1]-beta[m])*-1))
              }
             if (beta[m-1]!=0 & beta[m+1]==0){
                 lam2_1[m]=((-lambda2[m]*sign(beta[m]-beta[m-1]))- (lambda2[m]*sign(0-grad[m])*-1))
               }
                 }

         if(m==1){
           lam2_1[m]=-lambda2[m]*sign(beta[m+1]-beta[m])*-1
           }
          }


       if(m==length(beta)){
           if(beta[m-1]!=0){
           lam2_1[m]=-lambda2[m]*sign(beta[m]-beta[m-1])
              }
           if(beta[m-1]==0){
           lam2_1[m]=-lambda2[m]*sign(grad[m]-0)
           }
          }
        }


 if(length(i:m)!=1){
     if(nzb[i]) {

      cu_g1=cu_g2=cu_g
      cu_g1[i:m]=cumsum(grad[i:m])
      cu_g2[i:(m-1)]=cu_g1[m]-cu_g1[i:(m-1)]
   ##################################################
      P=numeric(length(cu_g1))
      P[i:m]=1
      P[i:m]=cumsum(P[i:m])
      P_grad1[i:m]=(cu_g1[i:m])/P[i:m]
      P_grad2[i:m]=(cu_g2[i:m])/(P[m]-P)[i:m]
      P_grad2[m]=0
    ##########################################

    ###############Calculate lambda1 and lambda2 with signs ################

     if(i==1){
      if (m!=length(beta)){
       lam2_1[i:m]=((-lambda2[i:m]*sign(P_grad2[i:m]-P_grad1[i:m])*-1))/P[i:m]
       lam2_1[m]=(-lambda2[m]*sign(beta[m+1]-beta[m])*-1)/(P[m])
       }
       if(m==length(beta)){
       lam2_1[i:m]=(-lambda2[i:m]*sign(P_grad2[i:m]-P_grad1[i:m])*-1)/P[i:m]
       lam2_1[m]=0
       }
       }

      if(i!=1){
       if(m!=length(beta)){
         lam2_1[i:m]=((-lambda2[i:m]*sign(beta[i]-beta[i-1])-(lambda2[i:m]*sign(P_grad2[i:m]-P_grad1[i:m])*-1)))/P[i:m]
         lam2_1[m]=(-lambda2[m]*sign(beta[i]-beta[i-1]) - (lambda2[m]*sign(beta[m+1]-beta[i])*-1))/(P[m])
         }
       if(m==length(beta)){
         lam2_1[i:m]=((-lambda2[i:m]*sign(beta[i]-beta[i-1])-(lambda2[i:m]*sign(P_grad2[i:m]-P_grad1[i:m])*-1)))/P[i:m]
         lam2_1[m]=(-lambda2[m]*sign(beta[i]-beta[i-1]))/(P[m])
         }
         }

     if(m!=length(beta)){
        lam2_2[i:m]=((-lambda2[i:m]*sign(P_grad2[i:m]-P_grad1[i:m])-(lambda2[i:m]*sign(beta[m+1]-beta[i])*-1)))/(P[m]-P[i:m])
        lam2_2[m]=0
        }
     if(m==length(beta)){
        lam2_2[i:m]=(-lambda2[i:m]*sign(P_grad2[i:m]-P_grad1[i:m]))/(P[m]-P[i:m])
        lam2_2[m]=0
      }
          #####################################################

      lam1_1[i:m] =(lambda1[i:m]*sign(P_grad1[i:m]))
      lam1_2[i:m]=(lambda1[i:m]*sign(P_grad2[i:m]))
      lam1_2[m]=0
      }


#######################################################################

  if(!nzb[i]){

      cu_g1=P_grad
      cu_g1[i:m]=cumsum(grad[i:m])
      cu_g2=numeric(length(cu_g1))
   ##################################################
      P=numeric(length(cu_g1))
      P[i:m]=1
      P[i:m]=cumsum(P[i:m])
      P_grad1[i:m]=(cu_g1[i:m])/P[i:m]
      P_grad2= numeric(length(beta))

    ##########################################

    ###############Calculate lambda1 and lambda2 with signs ################

     if(i==1){
      if (m!=length(beta)){
                         lam2_1[i:m]= (-lambda2[i:m]*sign(0-P_grad1[i:m])*-1)/P[i:m]
                         lam2_1[m]=(-lambda2[m]*sign(beta[m+1]-beta[m])*-1)/(P[m])
                        }
                       if(m==length(beta)){
                       lam2_1[i:m]= (-lambda2[i:m]*sign(0-P_grad1[i:m])*-1)/P[i:m]
                       lam2_1[m]=0
                       }
                       }

      if(i!=1){
             if (m == length(beta)){
                 if(beta[i-1]==0){
                  lam2_1[i:m]=(-lambda2[i:m]*sign(P_grad1[i:m]-0)- (lambda2[i:m]*sign(0-P_grad1[i:m])*-1))/P[i:m]
                  lam2_1[m]=(-lambda2[m]*sign(P_grad1[m]-0))/(P[m])
                  }
                  if(beta[i-1]!=0){
                  lam2_1[i:m]=((-lambda2[i:m]*(sign(beta[i:m]-beta[i-1])))- (lambda2[i:m]*sign(0-P_grad1[i:m])*-1))/P[i:m]
                  lam2_1[m]=(-lambda2[m]*(sign(beta[m]-beta[i-1])))/(P[m])
                  }
                 }
               if(m != length(beta)){
                 if(beta[m+1]==0 & beta[i-1]==0){
                   lam2_1[i:m]=((-lambda2[i:m]*(sign(P_grad1[i:m]-0)))- (lambda2[i:m]*sign(0-P_grad1[i:m])*-1))/P[i:m]
                   lam2_1[m]= lam2_1[m]/P[m]
                   }
                   if(beta[m+1]==0 & beta[i-1]!=0){
                   lam2_1[i:m]=((-lambda2[i:m]*(sign(beta[i:m]-beta[i-1])))- (lambda2[i:m]*sign(0-P_grad1[i:m])*-1))/P[i:m]
                   lam2_1[m]= lam2_1[m]/P[m]
                   }
                   if(beta[m+1]!=0 & beta[i-1]==0){
                   lam2_1[i:m]=((-lambda2[i:m]*(sign(P_grad1[i:m]-0)))- (lambda2[i:m]*sign(0-P_grad1[i:m])*-1))/P[i:m]
                   lam2_1[m]= ((-lambda2[m]*(sign(P_grad1[m]-0)))- (lambda2[m]*sign(beta[m+1]-beta[m])*-1))/P[m]

                   }
                   if(beta[m+1]!=0 & beta[i-1]!=0){
                   lam2_1[i:m]=((-lambda2[i:m]*(sign(beta[i:m]-beta[i-1])))- (lambda2[i:m]*sign(0-P_grad1[i:m])*-1))/P[i:m]
                   lam2_1[m] = ((-lambda2[m]*(sign(beta[m]-beta[i-1])))- (lambda2[m]*sign(beta[m+1]-beta[m])*-1))/P[m]

                   }
                }

              }
                         lam1_1[i:m]=(lambda1[i:m]*sign(P_grad1[i:m]))

        }


       }

     if(nzb[i] && (length(i:m)!=1)){

      direc1_1=numeric(length(cu_g1))
      direc1_2=numeric(length(cu_g2))

      direc1_1[i:m]=P_grad1[i:m]-lam1_1[i:m]+lam2_1[i:m]
      direc1_1[length(direc1_1)]=P_grad1[length(direc1_1)]-lam1_1[length(direc1_1)]+lam2_1[length(direc1_1)]
      direc1_2[i:m]=P_grad2[i:m]-lam1_2[i:m]+lam2_2[i:m]
      test=(direc1_1>direc1_2)
      test1=(direc1_1<direc1_2)
      testd = ((test&((P_grad1)>(P_grad2)))|(test1&((P_grad1)<(P_grad2))))
      testd[m] = F
      if(any(testd[i:m])){
      test2=(i:m)[which(testd[i:m])]

      direc_j=matrix(0,nrow=length(direc1_1),ncol=length(direc1_2))
      if(length(test2)>1){
      for(ik in 1:length(test2)){direc_j[i:test2[ik],test2[ik]]=(direc1_1)[test2[ik]]}
      for(ik in 1:(length(test2)-1)){direc_j[(test2[ik]+1):m,test2[ik]]=(direc1_2)[test2[ik]]}
      if(test2[length(test2)]!=m){
      direc_j[(test2[length(test2)]+1):m,test2[length(test2)]]=direc1_2[test2[length(test2)]]
      }
      } else
      if(length(test2)==1){
      if (test2==m){
      direc_j[i:test2,test2]=(direc1_1)[test2]
      }
      if (test2!=m){
      direc_j[i:test2,test2]=(direc1_1)[test2]
      direc_j[(test2+1):m,test2]=(direc1_2)[test2]
      }
      }
      norm1 <- apply(direc_j,2,function(il){sqrt(sum(t(il)*il))})
      vmax2=max(norm1)


      direc_i[i:m]=direc_j[(i:m),which(norm1==vmax2)]
      if(which(vmax2==norm1)==m){vmax=0}
      if(which(vmax2==norm1)!=m){
      vmax= sqrt(sum(t(direc_i)*direc_i)) }
      }
      }
  if(!nzb[i]){
              if(length(i:m)!=1){
              lam=numeric(length(lam1_1))
               lam=((lam1_1)-(lam2_1))

                if(any((P_grad1[i:m]*sign(P_grad1[i:m]))>(lam[i:m]*sign(P_grad1[i:m])))){
                 d_i = which((P_grad1[i:m]*sign(P_grad1[i:m]))>(lam[i:m]*sign(P_grad1[i:m])))
                 direc1_1= P_grad1 - lam1_1 + lam2_1

                 direc_j=matrix(0,length(i:m),length(i:m))
                 for(im in 1:length(d_i)){
                 direc_j[1:d_i[im],d_i[im]]=((P_grad1[i:m])[d_i[im]]-((lam1_1[i:m])[d_i[im]])+((lam2_1[i:m])[d_i[im]]))
                 }
                 norm_d <- apply(direc_j,2,function(il){sqrt(sum(t(il)*il))})
                 vmax1=max(norm_d)
                 direc1_1[i:m]=direc_j[,which(norm_d==vmax1)]
                 direc_i[i:m]=direc1_1[i:m]

                 vmax = sqrt(sum(t(direc_i)*direc_i))

                 }

              }
             if(length(i:m)==1){

                   if((grad[m]*sign(grad[m]))>(lam1_1[m]-lam2_1[m])*sign(grad[m])){
                      direc_i[m]=grad[m] - lam1_1[m] + lam2_1[m]
                       vmax= sqrt(sum(t(direc_i)*direc_i))
                      }
                   }
              }

                           if (oldvmax>=vmax) {direc_i=olddirec_i
                               vmax=oldvmax}


                      olddirec_i=direc_i
                      oldvmax=vmax

          if(nzb[i]){
                     if(length(i:m)==1){io=i+1}
                            if(length(i:m)!=1){
                                if(j==m){io=m+1}
                                  if(j!=m){io=j+1}
                                              }
                                       i=io  } else {i=i+1}


         }
}






    return(list(direc_i = direc_i))

}

######################################

######################################


 tryNRstep <- function (beta = beta,direc=direc,active=active,X=X,localfit=localfit,newfit=newfit)   {

 NRm=nonZero(beta)
                 if(any((beta==0) & active)){
                    NRact=which((beta==0) & active)
                    NRm=rbind(NRm,cbind(NRact,NRact))
                    }
                 NRM=rbind(NRm,c(0,0))
                 NR_X=matrix(0,nrow=nrow(X),ncol=sum(active))
                 NRb=beta[active]
                 NR_X = as.matrix(X[,active])


                 if (is.list(localfit$W)) {
                   XdW <- NR_X * matrix(sqrt(localfit$W$diagW), nrow(NR_X), ncol(NR_X))
                   hessian <- -crossprod(XdW) + crossprod(crossprod(localfit$W$P, NR_X))
                   } else if (length(localfit$W) > 1) {
                   XW <- NR_X * matrix(sqrt(localfit$W), nrow(NR_X), ncol(NR_X))
                   hessian <- -crossprod(XW)
                   } else {
                   hessian <- -crossprod(NR_X)
                   }
                   NRbeta <- NRb - drop(.solve(hessian, direc[active]))
                   newbeta=numeric(length(beta))
                   newbeta[active] = NRbeta
                   NRfailed <- !all(sign(newbeta[beta!=0]) == sign(beta[beta!=0]))
          if (!NRfailed) {
            NR_X=matrix(0,nrow=nrow(X),ncol=(nrow(NRM)-1))
                 NRb=beta[NRM[-(nrow(NRM)),1]]
                 for (ik in 1:(nrow(NRM)-1)){if(NRM[ik,1]==NRM[ik,2]){
                   NR_X[,ik]=X[,NRM[ik,1]]}
                   else{NR_X[,ik]=rowSums(X[,NRM[ik,1]:NRM[ik,2]])}
                  }
                 NR_X=as.matrix(NR_X)
          if (is.list(localfit$W)) {
                   XdW <- NR_X * matrix(sqrt(localfit$W$diagW), nrow(NR_X), ncol(NR_X))
                   hessian <- -crossprod(XdW) + crossprod(crossprod(localfit$W$P, NR_X))
                   } else if (length(localfit$W) > 1) {
                   XW <- NR_X * matrix(sqrt(localfit$W), nrow(NR_X), ncol(NR_X))
                   hessian <- -crossprod(XW)
                   } else {
                   hessian <- -crossprod(NR_X)
                   }
                   NRbeta <- NRb - drop(.solve(hessian, direc[NRM[-(nrow(NRM)),1]]))
                   newbeta=numeric(length(beta))
                   for (ki in 1:(nrow(NRM)-1)){newbeta[NRM[ki,1]:NRM[ki,2]]=NRbeta[ki]}
                   beta=newbeta
          newfit <- TRUE
        }
        return(list(beta=beta,newfit=newfit,NRfailed = NRfailed))

        }

##############################################

