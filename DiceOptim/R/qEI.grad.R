decompo <- function(x, mu, Sigma,k) {
  ##### 
  if (length(k)==length(x[,1])) {
      return(dmnorm(c(x),c(mu),Sigma))
  }
  
  x1 <- x[k,1]
  x2 <- x[-k,1]
  mu1 <- mu[k,1]
  Sig22 <- Sigma[-k,-k]
  Sig21 <- matrix(Sigma[-k,k],ncol=length(k))
  Sig12 <- t(Sig21)
  Sig11 <- Sigma[k,k]
  Sig11Inv <- solve(Sig11)
  varcov <- Sig22 - Sig21%*%Sig11Inv%*%Sig12
  varcov <- 0.5*(varcov+t(varcov))
  if (min(diag(varcov))<=0) {diag(varcov) <- diag(varcov)*(diag(varcov)>=0)+10^(-9)}
  moy <- mu[-k,1]+Sig21 %*% Sig11Inv%*%(x1-mu1)

  
  return (
    dmnorm(x=x1,  mu[k,1], varcov = Sig11)
    *
    pmnorm(x=t(matrix(x2)), mean = t(moy), varcov = varcov,maxpts=(length(moy)*200))
  )
}

# gradient of CDF
GPhi <- function(x,mu,Sigma) {
  q <- length(x[,1])
  res <- matrix(NaN,q,1)
  for (i in 1:q) {
    res[i,1] <- decompo(x,mu,Sigma,i)
  }
  return( res )
}

# hessian matrix of CDF
HPhi <- function(x,mu,Sigma,gradient=NULL) {
  q <- length(x[,1])
  if (q == 1) {
    res <- -(x-mu)/Sigma*dnorm(x,mu,sqrt(Sigma))
  }
  else {
    res <- matrix(0,q,q)
    for (i in 1:(q-1)) { 
      for (j in (i+1):q) {
        res[i,j] <- decompo(x,mu,Sigma,c(i,j))
      }
    }
    res <- res + t(res) # hessian matrix is symetric
    # diagonal terms can be computed with the gradient of CDF and the other hessian terms
    if(is.null(gradient)) {
      res <- res - diag(((c(x)-c(mu))*c(GPhi(x,mu,Sigma)) + diag(Sigma%*%res))/diag(Sigma))
    } else {
      res <- res - diag(((c(x)-c(mu))*c(gradient) + diag(Sigma%*%res))/diag(Sigma))
    }
  }
  return(res)
}

krigingDeriv <- function(x, model, type="UK", envir=NULL){

# Compute from a km object the derivatives of kriging mean and covariance in respect 
# of a batch of points.

  ########################################################################################
  q <- length(x[,1])
  ########################################################################################
  # Convert x in proper format(s)
  d <- length(x[1,])
  if (d != model@d){ stop("x does not have the right size") }
  newdata.num <- matrix(as.numeric(x),ncol=d)
  newdata <- data.frame((newdata.num))
  colnames(newdata) = colnames(model@X)
  ########################################################################################
  # Get quantities related to the model
  T <- model@T
  X <- model@X
  z <- model@z
  u <- model@M
  covStruct <- model@covariance
  F.newdata <- model.matrix(model@trend.formula, data=newdata)
  
  # Get quantities related to the prediction
  if (is.null(envir))
  {  
    predx <- predict(object=model, newdata=newdata, type=type,
                     checkNames = FALSE,se.compute=FALSE,cov.compute=TRUE)
    kriging.mean <- predx$mean
    kriging.cov <- predx$cov
    v <- predx$Tinv.c
    c <- predx$c
    
  } else
  {  # If uploaded through "envir", no prediction computation is necessary 
    toget <- matrix(c("kriging.mean","kriging.cov","Tinv.c","c"), 4, 1)
    apply(toget, 2, get, envir=envir)
    kriging.mean <- envir$kriging.mean
    kriging.cov  <- envir$kriging.cov
    v            <- envir$Tinv.c
    c            <- envir$c
  }

  ########################################################################################
  # Pursue calculation only if standard deviation is non-zero
  if ( prod(diag(kriging.cov)/(model@covariance@sd2) < 1e-08) ) 
  { return (NULL)
  } else 
  { 
    kriging.mean.jacob <- array(0,dim=c(q,d,q))  # three dimensions : d(m_k)/d(x_ij)
    kriging.cov.jacob <- array(0,dim=c(q,d,q,q)) # four dimensions : d(Sigma_kl)/d(x_ij)
    tuuinv <- solve(t(u)%*%u)
    covM <- covMatrix(covStruct, as.matrix(newdata))[[1]]
    
    for(k in 1:q) {
    # Compute derivatives of the covariance and trend functions
    dc <- covVector.dx(x=newdata.num[k,], X=X, object=covStruct, c=c[,k])  
    f.deltax <- trend.deltax(x=newdata.num[k,], model=model)
    # Compute gradients of the kriging mean at point k
    W <- backsolve(t(T), dc, upper.tri=FALSE)
    kriging.mean.jacob[k,,k] <- t(z)%*%W + model@trend.coef%*%f.deltax
    # Compute gradients of the kriging covariance between point k and point l
    for (l in 1:q) {
        # Compute gradients of the kernel between point k and point l
        ker.grad <- covVector.dx(x=newdata.num[k,], X=newdata[l,],
                                 object=covStruct, c=covM[k,l])
        kriging.cov.jacob[k,,k,l] <- ker.grad - t(v[,l])%*%W
        if (type=="UK") {
          kriging.cov.jacob[k,,k,l] <- kriging.cov.jacob[k,,k,l] + t((t(f.deltax) - t(W)%*%u)%*%tuuinv%*%t(F.newdata[l,] - t(v[,l])%*%u))
        }
      }
    }
  }
  for(l in 1:q) {
    for(j in 1:d) {
       kriging.cov.jacob[l,j,,] <- kriging.cov.jacob[l,j,,] + t(kriging.cov.jacob[l,j,,])
    }
  }
  ########################################################################################  
  return(list(kriging.mean,kriging.mean.jacob,kriging.cov,kriging.cov.jacob))
}


qEI.grad <- function(x, model, plugin=NULL, type="UK", minimization = TRUE, fastCompute = TRUE, eps = 10^(-6), envir=NULL) {
  
  
  #compteur_global <<- compteur_global + 1
  
  if (!minimization) {
    if (is.null(plugin)) {
      plugin <- -max(model@y)
    } 
  }
  
  # Compute the multipoint expected improvement from a batch of point and a km object.
  #### Arguments:
  # x: position of multEI evaluation. One row correspond to one point.
  # model: a km object.
  # plugin: the threshold used to compute the EI
  # type: do we perform a universal krigiging ("UK") ?
  d <- model@d
  if (!is.matrix(x)) {
    x <- matrix(x,ncol = d)
  }

  q <- nrow(x)
  xb <- rbind(model@X,x)
  ux <- unique(round(xb,digits=8))
  if(length(xb[,1])!=length(ux[,1])) {return (matrix(0,q,d))}
  
  if (is.null(plugin)) plugin <- min(model@y)


  if (q == 1) {
    if (!minimization) {
      stop("qEI.grad doesn't work in \'minimization = FALSE\' when dim = 1 (in progress).")
    }
    return(EI.grad(x,model,plugin,type,envir))
  }
  
  if(!is.null(envir)) {
    
    if (fastCompute == TRUE) {
      toget <- matrix(c("pk"), 1, 1)
      apply(toget, 2, get, envir=envir)
      pk <- envir$pk
    } else {
      toget <- matrix(c("pk","symetric_term"), 2, 1)
      apply(toget, 2, get, envir=envir)
      pk <- envir$pk
      symetric_term <- envir$symetric_term
   }
  }
  
  # Get kriging (conditionnal) mean and covariance and their derivatives
  kriging.dx <- krigingDeriv(x=x, model=model, type=type, envir=envir)
  if (is.null(kriging.dx)) {return (matrix(0,q,d))}
  kriging.mean        <- kriging.dx[[1]] # 1 dimension
  kriging.mean.jacob  <- kriging.dx[[2]] # 3 dimensions q*d*q,  kriging.mean.jacob[l,j,k] = d(m_k)/d(x_lj)
  if (!minimization) {kriging.mean <- -kriging.mean; kriging.mean.jacob <- -kriging.mean.jacob}
  kriging.cov         <- kriging.dx[[3]] # 2 dimensions q*q
  kriging.cov.jacob   <- kriging.dx[[4]] # 4 dimensions q*d*q*q : kriging.cov.jacob[l,j,k,i]= d(Sigma_ki)/d(x_lj)
  # Initialisation
  EI.grad <- matrix(0,q,d) # result
  b <- matrix(rep(0,q),q,1)
  L <- -diag(rep(1,q))
  Dpk <- matrix(0,q,d)
  if (fastCompute == TRUE) {
    termB <- matrix(0,q,d)
  }
  Sigk_dx <- array(0,dim=c(q,d,q,q))
  mk_dx <- array(0,dim=c(q,d,q))
  # First sum of the formula
  for (k in 1:q) {
    bk <- b
    bk[k,1] <- plugin     # creation of vector b^(k)
    Lk <- L
    Lk[,k] <- rep(1,q)    # linear application to transform Y to Zk
    tLk <- t(Lk)
    mk     <- Lk %*% kriging.mean # mean of Zk (written m^(k) in the formula)
    Sigk <- Lk %*% kriging.cov %*% tLk # covariance of Zk
    Sigk <- 0.5*(Sigk+t(Sigk)) # numerical symetrization
    
    # term A1
    if (is.null(envir)) {
      EI.grad[k,] <-  EI.grad[k,] - kriging.mean.jacob[k,,k]*pmnorm(x=t(bk), mean = t(mk), varcov = Sigk,maxpts=q*200)
      } else {
      EI.grad[k,] <-  EI.grad[k,] - kriging.mean.jacob[k,,k]*pk[k]
    }
    
    # compute gradient ans hessian matrix of the CDF term pk.
    gradpk <- GPhi(bk-mk,b,Sigk)
    hesspk <- HPhi(x = bk, mu = mk, Sigma = Sigk, gradient = gradpk)

    # term A2
    for (l in 1:q) {
      for (j in 1:d) {
        Sigk_dx[l,j,,] <- Lk%*%kriging.cov.jacob[l,j,,]%*%tLk
        mk_dx[l,j,] <- Lk%*%kriging.mean.jacob[l,j,]
        Dpk[l,j] <- 0.5*sum(hesspk*Sigk_dx[l,j,,]) - crossprod(gradpk,mk_dx[l,j,])
      }
    }
    EI.grad <- EI.grad + (plugin-kriging.mean[k])*Dpk
    # term B
    if (fastCompute == TRUE) {
      gradpk1<- GPhi(bk-mk+Sigk[,k]*eps,b,Sigk)
      hesspk1<- HPhi(x = bk+Sigk[,k]*eps, mu = mk, Sigma = Sigk, gradient = gradpk1)
      for (l in 1:q) {
	for (j in 1:d) {
	  f1 <- -crossprod(mk_dx[l,j,],gradpk1) + eps*crossprod(Sigk_dx[l,j,,k],gradpk1) + 0.5*sum(Sigk_dx[l,j,,]*hesspk1)
	  f  <- -crossprod(mk_dx[l,j,],gradpk ) + 0.5*sum(Sigk_dx[l,j,,]*hesspk)
	  termB[l,j] <- 1/eps*(f1-f)
	}
      }
    } else {
      B1 <- B2 <- B3 <- matrix(0,q,d)
      for (i in 1:q) {
	Sigk_ik <- Sigk[i,k]
	Sigk_ii <- Sigk[i,i]
	mk_i <- mk[i]
	mk_dx_i <- mk_dx[,,i]
	bk_i <- bk[i,1]
	ck_pi <- bk[-i,1]-mk[-i,]-(bk[i,1]-mk[i,1])/Sigk_ii*Sigk[-i,i]
	Sigk_pi <- 0.5*(Sigk[-i,-i] - 1/Sigk_ii*Sigk[-i,i]%*%t(Sigk[-i,i])+t(Sigk[-i,-i] - 1/Sigk_ii*Sigk[-i,i]%*%t(Sigk[-i,i])))
	Sigk_dx_ii <- Sigk_dx[,,i,i]
	Sigk_dx_ik <- Sigk_dx[,,i,k]
	phi_ik <- dnorm(bk[i,1],mk_i,sqrt(Sigk_ii))
	dphi_ik_dSig <- ((bk_i-mk_i)^2/(2*Sigk_ii^2)-0.5*1/Sigk_ii)*phi_ik
	dphi_ik_dm <- (bk_i-mk_i)/Sigk_ii*phi_ik
	if (is.null(envir)) {
	  Phi_ik <- pmnorm(x=ck_pi, mean = rep(0,q-1), varcov = Sigk_pi,maxpts=(q-1)*200)
	} else {
	  Phi_ik <- symetric_term[k,i]/phi_ik
	}
	GPhi_ik <- GPhi(matrix(ck_pi,q-1,1),matrix(0,q-1,1),Sigk_pi)
	HPhi_ik <- HPhi(matrix(ck_pi,q-1,1),matrix(0,q-1,1),Sigk_pi,GPhi_ik)
	Sigk_mi <- Sigk[-i,i]
	for (l in 1:q) {
	  for (j in 1:d) {
	    # B1
	    B1[l,j] <- B1[l,j] + Sigk_dx_ik[l,j]*phi_ik*Phi_ik
	    # B2
	    B2[l,j] <- B2[l,j] + Sigk_ik*(mk_dx_i[l,j]*dphi_ik_dm + dphi_ik_dSig*Sigk_dx_ii[l,j])*Phi_ik
	    # B3
	    dck_pi <- -mk_dx[l,j,-i]+(mk_dx_i[l,j]*Sigk_ii+(bk_i-mk_i)*Sigk_dx_ii[l,j])/Sigk_ii^2*Sigk_mi-(bk[i,1]-mk_i)/Sigk_ii*Sigk_dx[l,j,-i,i]
	    SigtCross <- tcrossprod(Sigk_dx[l,j,-i,i],Sigk_mi)
	    dSigk_pi <- Sigk_dx[l,j,-i,-i]+Sigk_dx_ii[l,j]/Sigk_ii^2*tcrossprod(Sigk_mi,Sigk_mi)-Sigk_ii^-1*(SigtCross+t(SigtCross))
	    B3[l,j] <- B3[l,j] + Sigk_ik*phi_ik*(crossprod(GPhi_ik,dck_pi)+0.5*sum(HPhi_ik*dSigk_pi))
	  }
	}
      }
    }
    if (fastCompute == TRUE) {
      EI.grad <- EI.grad + termB
    } else {
      EI.grad <- EI.grad + B1 + B2 + B3
    }
  }
  if(is.nan(sum(EI.grad))) {return(matrix(0,q,d))}
  return (EI.grad)
}
