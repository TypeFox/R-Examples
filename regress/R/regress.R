regress <- function(formula, Vformula, identity=TRUE, kernel=NULL,
                    start=NULL, taper=NULL, pos, verbose=0, gamVals=NULL,
                    maxcyc=50, tol=1e-4, data){

  ## Vformula can just be something like ~ V0 + V1
  ## or leave it out or Vformula=NULL
  ## assume its in the form ~ V1 + V2 + ... + Vn or missing or Vformula=NULL
  ## for random effects and random interactions for factors A and B include
  ## ~ A + B + I(A:B)

  if(verbose>9) cat("Extracting objects from call\n")
  if(missing(data)) data <- environment(formula)
  mf <- model.frame(formula,data=data,na.action=na.pass)
  mf <- eval(mf,parent.frame())
  y <- model.response(mf)

  model <- list()
  model <- c(model,mf)

  if(missing(Vformula)) Vformula <- NULL

  ## Find missing values in fixed part of the model :: Change Aui Mar 1 2012
  isNA <-  apply(is.na(mf), 1, any)

  if(!is.null(Vformula))
  {
      V <- model.frame(Vformula,data=data,na.action=na.pass)
      V <- eval(V, parent.frame())
      # find missings in random part of the model Aui Mar 1 2012
      mfr <- is.na(V)
      if(ncol(mfr) == 1){
        isNA <- isNA | mfr
      } else {
        isNA <- isNA | apply(mfr[,!apply(mfr,2,all)], 1, any) # use only columns of the matrix with some nonmissing values
      }
      rm(mfr)
      Vcoef.names <- names(V)
      V <- as.list(V)
      k <- length(V)
  } else {
      V <- NULL
      k <- 0
      Vcoef.names=NULL
  }

  if(ncol(mf)==1) mf <- cbind(mf,1)
  X <- model.matrix(formula, mf[!isNA,]) # Aui Mar 1 2012 account for missings in the random part

  y <- y[!isNA]
  n <- length(y)
  Xcolnames <- dimnames(X)[[2]]
  if(is.null(Xcolnames)) {
      Xcolnames <- paste("X.column",c(1:dim(as.matrix(X))[2]),sep="")
  }

  X <- matrix(X, n, length(X)/n)
  qr <- qr(X)
  rankQ <- n-qr$rank
  if(qr$rank) {
      X <- matrix(X[, qr$pivot[1:qr$rank]],n,qr$rank)
      Xcolnames <- Xcolnames[qr$pivot[1:qr$rank]]
  } else {
      cat("\nERROR: X has rank 0\n\n")
  }

  if(verbose>9) cat("Setting up kernel\n")

  if(missing(kernel)){
      K <- X
      colnames(K) <- Xcolnames
      reml <- TRUE
      kernel <- NULL
  } else {
      if(length(kernel)==1 && kernel>0){
          K <- matrix(rep(1, n), n, 1)
          colnames(K) <- c("1")
      }
      if(length(kernel)==1 && kernel<=0) {
          K <- Kcolnames <- NULL
          KX <- X
          rankQK <- n
      }
      if(length(kernel) > 1) {
          ##K is a matrix I hope :: Change Aui Mar 1 2012
          if(is.matrix(kernel)) {
              K <- kernel[!isNA,]
          } else {
              K <- model.frame(kernel, data=data, na.action=na.pass)
              K <- eval(K, parent.frame())
              if(ncol(K) == 1){
                dimNamesK <- dimnames(K)
                K <- K[!isNA, ]
                dimNamesK[[1]] <- dimNamesK[[1]][!isNA]
                K <- data.frame(V1 = K)
                dimnames(K) <- dimNamesK
              } else {
                K <- K[!isNA, ]
              }
              K <- model.matrix(kernel, K)
          }
      }
      reml <- FALSE
  }

  if(!is.null(K)){
      Kcolnames <- colnames(K)
      qr <- qr(K)
      rankQK <- n - qr$rank
      if(qr$rank == 0) K <- NULL else {
          K <- matrix(K[, qr$pivot[1:qr$rank]],n,qr$rank)
          Kcolnames <- Kcolnames[qr$pivot[1:qr$rank]]
          KX <- cbind(K, X) # Spanning K + X: Oct 12 2011
          qr <- qr(KX)
          KX <- matrix(KX[, qr$pivot[1:qr$rank]],n,qr$rank) # basis of K+X
      }
  }

  if(missing(maxcyc)) maxcyc <- 50
  if(missing(tol)) tol <- 1e-4
  delta <- 1

  if(verbose>9) cat("Removing parts of random effects corresponding to missing values\n")
  ## remove missing values
  for(i in 1:k)
  {
      if(is.matrix(V[[i]]))
      {
          V[[i]] <- V[[i]][!isNA, !isNA]
      }
      if(is.factor(V[[i]]))
      {
          V[[i]] <- V[[i]][!isNA]
      }
  }

  In <- diag(rep(1,n),n,n)

  if(identity) {
      V[[k+1]] <- as.factor(1:n)
      names(V)[k+1] <- "In"
      k <- k+1

      Vcoef.names <- c(Vcoef.names,"In")
      Vformula <- as.character(Vformula)
      Vformula[1] <- "~"
      Vformula[2] <- paste(Vformula[2],"+In")
      Vformula <- as.formula(Vformula)
  }

  model <- c(model,V)
  model$formula <- formula
  model$Vformula <- Vformula

  ## specify which parameters are positive and which are negative
  ## pos = c(1,1,0) means first two parameters are positive, third is either
  if(!missing(pos)) pos <- as.logical(pos)
  if(missing(pos)) pos <- rep(FALSE,k)
  if(length(pos) < k) cat("Warning: argument pos is only partially specified; additional terms (n=",k-length(pos),") set to FALSE internally.\n",sep="")
  pos <- c(pos,rep(FALSE,k))
  pos <- pos[1:k]

  ## Sherman Morrison Woodbury identities for matrix inverses can be brought to bear here
  if(verbose>9) cat("Checking if we can apply the Sherman Morrison Woodbury identites for matrix inversion\n")
  if (all(sapply(V, is.factor)) & k>2 ) {  # Contribution by Hans Jurgen Auinger
      SWsolveINDICATOR <- TRUE
  } else SWsolveINDICATOR <- FALSE
  Z <- list()
  for (i in 1:length(V)) {
      if (is.factor(V[[i]])) {
          Vi <- model.matrix(~V[[i]] - 1)
          colnames(Vi) <- levels(V[[i]])
          Z[[i]] <- Vi
          V[[i]] <- tcrossprod(Vi)
      } else{
          Z[[i]] <- V[[i]]
      }
  }
  names(Z) <- names(V)

  ## So V is always a list of variance coavriance matrices, Z contains
  ## the model matrices of factors when we need to invoke the Sherman
  ## Woodbury identities

  ## Expected Fisher Information
  A <- matrix(rep(0, k^2), k, k)
  entries <- expand.grid(1:k,1:k)

  x <- rep(0,k)
  sigma <- c(1,rep(0, k-1))
  stats <- rep(0, 0)

  ## START ALGORITHM

  if(missing(taper)){
      taper <- rep(0.9, maxcyc)
      if(missing(start) && k>1) taper[1:2] <- c(0.5, 0.7)
  } else {
      taper <- pmin(abs(taper), 1)
      if((l <- length(taper)) < maxcyc) taper <- c(taper, rep(taper[l], maxcyc-l))
  }

  if(!is.null(start)) {
      ## pad start with zeros if required
      start <- c(start, rep(1,k))
      start <- start[1:k]
  }

  if(k>2 && is.null(start)) start <- rep(var(y,na.rm=TRUE),k)
  if(k==1 && is.null(start)) start <- var(y,na.rm=TRUE)

  if(is.null(start) && k==2) {
      if(missing(gamVals)) {
          gamVals <- seq(0.01,0.02,length=3)^2
          gamVals <- sort(c(gamVals,seq(0.1,0.9,length=3),1-gamVals))
          gamVals <- 0.5
      }
      if(length(gamVals)>1) {
          if(verbose>=1) cat("Evaluating the llik at gamma = \n")
          if(verbose>=1) cat(gamVals)
          if(verbose>=1) cat("\n")
          reg.obj <- reml(gamVals,y,X,V[[1]],V[[2]],verbose=verbose)
          llik <- reg.obj$llik
          llik <- as.double(llik)
          if(verbose>=2) cat(llik,"\n")
          gam <- gamVals[llik==max(llik)]
          gam <- gam[1]
          if(verbose>=2) cat("MLE is near",gam,"and llik =",max(llik),"there\n")
      }
      if(length(gamVals)==1) {
          ## go straight to the Newton Raphson at gamVals
          gam <- gamVals[1]
          reg.obj <- list(rms=var(y))
      }
      start <- c(1-gam,gam)*reg.obj$rms
      ## it tends to take huge steps when starting at gam=0.9999
      if(gam==0.9999) {
          taper[1] <- taper[1]/100
          maxcyc <- maxcyc*10
      }
      if(verbose>=1) cat(c("start algorithm at",round(start,4),"\n"))
  }

  if(is.null(start) & k>2) {
      ## Never gets here by default - but this could be implemented,
      ## though it does add on a few extra iterations at the
      ## start.... not necessary in basic examples

      LLvals <- NULL
      ## equal weights
      V2 <- V[[2]]
      for(ii in 3:k) V2 <- V2 + V[[ii]]
      LLvals <- c(LLvals,reml(0.5,y,X,V[[1]],V2)$llik)
      ## Most at one end
      V2 <- V[[1]] + V2 ## total of all Vs
      for(ii in 1:k) {
          V2 <- V2 - V[[ii]]
          LLvals <- c(LLvals,reml(0.75,y,X,V2,V[[ii]])$llik)
      }
      best <- which.max(LLvals)
      if(verbose) {
          cat("Checking starting points\n")
          cat("llik values of", LLvals, "\n")
      }
      if(best==1) {
          start <- rep(var(y,na.rm=TRUE),k)
      } else {
          start <- rep(0.25,k)
          start[best] <- 0.75
      }
  }

  sigma <- coef <- start
  ## reparameterise so everything will get into the correct spot after exp
  coef[pos] <- log(sigma[pos])
  coef[!pos] <- sigma[!pos]

  ## Set the memory requirements beforehand
  T <- vector("list", length=k)
  for(ii in 1:k) T[[ii]] <- matrix(NA,n,n)

  for(cycle in 1:maxcyc){

      ## Limit how far out we go on the logarithmic scale
      ind <- which(pos)
      if(length(ind)) {
          coef[ind] <- pmin(coef[ind],20)
          coef[ind] <- pmax(coef[ind],-20) ## so on regular scale everything is between exp(-20) and exp(20)
          sigma[ind] <- exp(coef[ind])
      }

      if(verbose) {
          cat(cycle, "sigma =",sigma)
          ##cat(sigma)
      }

      if(!SWsolveINDICATOR) {
          Sigma <- 0
          ## can we get rid of this loop?
          for(i in 1:k) Sigma <- Sigma + V[[i]]*sigma[i]

          W <- solve(Sigma,In)
      } else {
          W <- SWsolve2(Z[1:(k-1)],sigma)
      }

      if(is.null(K)) WQK <- W else {
          WK <- W %*% K
          WQK <- W - WK %*% solve(t(K)%*%WK, t(WK))
      }
      if(reml) WQX <- WQK else {
          WX <- W %*% KX        # including the kernel (Oct 12 2011)
          WQX <- W - WX %*% solve(t(KX)%*%WX, t(WX))
      }

      rss <- as.numeric(t(y) %*% WQX %*% y)

      ##if(verbose>9) cat("Sigma[1:5]",Sigma[1:5],"\n")
      ##if(verbose>9) cat("RSS",rss,"WQX[1:5]",WQX[1:5],"\n")
      sigma <- sigma * rss/rankQK
      coef[!pos] <- sigma[!pos]
      coef[pos] <- log(sigma[pos])
      WQK <- WQK * rankQK/rss
      WQX <- WQX * rankQK/rss
      rss <- rankQK ## looks bad but the rss is absorbed into WQK so the rss term comes out of eig below

      eig <- sort(eigen(WQK,symmetric=TRUE,only.values=TRUE)$values, decreasing=TRUE)[1:rankQK]
      if(any(eig < 0)){
          cat("error: Sigma is not positive definite on contrasts: range(eig)=", range(eig), "\n")
          WQK <- WQK + (tol - min(eig))*diag(nobs)
          eig <- eig + tol - min(eig)
      }
      ldet <- sum(log(eig))
      llik <- ldet/2 - rss/2
      if(cycle == 1) llik0 <- llik
      delta.llik <- llik - llik0
      llik0 <- llik

      ## From Jean-Luc Jannick to fix excess carriage returns
      if (verbose){
        if (reml) cat(" resid llik =", llik, "\n") else cat(" llik =", llik, "\n")
        cat(cycle, "adjusted sigma =", sigma)
        if (cycle > 1){
          if (reml) cat(" delta.llik =", delta.llik, "\n") else cat(" delta.llik =", delta.llik, "\n")
        } else cat("\n")
      }
      
      ## now the fun starts, derivative and expected fisher info
      ## the 0.5 multiple is ignored, it is in both and they cancel

      ##T <- list(NULL)
      x <- NULL

      ## derivatives are now D[[i]] = var.components[i]*V[[i]]
      var.components <- rep(1,k)
      ind <- which(pos)
      if(length(ind)) var.components[ind] <- sigma[ind]

      ## Slow part - order k n-squared
      if(!SWsolveINDICATOR) {
          if(identity) {
              T[[k]] <- WQK
              if(k>1) {
                  for(ii in (k-1):1) T[[ii]] <- WQK %*% V[[ii]]
              }
          } else {
              for(ii in 1:k) T[[ii]] <- WQK %*% V[[ii]]
          }
      } else {
          if(identity) {
              T[[k]] <- WQK
              if(k>1) {
                  for(ii in (k-1):1) T[[ii]] <- tcrossprod(WQK %*% Z[[ii]],Z[[ii]])
              }
          } else {
              for(ii in 1:k) T[[ii]] <- tcrossprod(WQK %*% Z[[ii]],Z[[ii]])
          }
      }


      x <- sapply(T,function(x) as.numeric(t(y) %*% x %*% WQX %*% y - sum(diag(x))))
      x <- x * var.components

      ## See nested for loops commented out below - evaluating the Expected Fisher Information, A
      ff <- function(x) sum(T[[x[1]]] * t(T[[x[2]]])) * var.components[x[1]] * var.components[x[2]]
      aa <- apply(entries,1,ff)
      A[as.matrix(entries)] <- aa

      ##for(i in 1:k)
      ##{
      ##if(identity && i==1 && pos[i]) {
      ##  T[[i]] <- WQ
      ##} else

      ##if(SWsolveINDICATOR==FALSE) {
      ##    T[[i]] <- WQ %*% V[[i]]
      ##} else {
      ##    T[[i]] <- WQ %*% tcrossprod(V[[i]])
      ##}
      ##x[i] <- as.numeric(t(y) %*% T[[i]] %*% WQ %*% y - sum(diag(T[[i]])))
      ##x[i] <- x[i]*var.components[i]

      ##for(j in c(1:i))
      ##  {
      ## expected fisher information
      ## ##A[j,i] <- Tr(T[[j]] %*% T[[i]])
      ## ##the Ts are not symmetric, hence the transpose below

      ##A[j,i] <- sum(T[[j]] * t(T[[i]]))
      ##A[j,i] <- A[j,i]*var.components[i]*var.components[j]
      ##A[i,j] <- A[j,i]
      ##}
      ##}

      stats <- c(stats, llik, sigma[1:k], x[1:k])
      if(verbose>=9) {
          ##cat(c(rllik1, rllik2, sigma[1:k], x[1:k]),"\n")
      }

      A.svd <- ginv(A)
      x <- A.svd %*% x

      if(qr(A)$rank < k){
          if(cycle==1) {
              if(verbose) {
                  cat("Warning: Non identifiable dispersion model\n")
                  ##print(round(A,6))
                  cat(sigma)
                  cat("\n")
              }
          }
      }

      ## end of newton-raphson step
      ## x is  -l'(sigma)/l''(sigma)
      ## hence we add instead of subtract
      ## taper controls the proportion of each step to take
      ## for some reason 0.5 works very well

      ##if(all(pos==1)) x <- sign(x) * pmin(abs(x),5) ## limit maximum shift we can take in one step to 5 units on log scale
      coef <- coef + taper[cycle] * x
      sigma[!pos] <- coef[!pos]
      sigma[pos] <- exp(coef[pos])

      ## check the change in llik is small
      if(cycle > 1 & abs(delta.llik) < tol*10) break
      if(max(abs(x)) < tol) break
  }

  ## Recompute Sigma at adjusted sigman values
  if(!SWsolveINDICATOR) {
      Sigma <- 0
      ## can we get rid of this loop?
      for(i in 1:k) Sigma <- Sigma + V[[i]]*sigma[i]

      W <- solve(Sigma,In)
  } else {
      W <- SWsolve2(Z[1:(k-1)],sigma)
  }


  if(cycle==maxcyc)
  {
      ## issue a warning
      if(verbose) cat("WARNING:  maximum number of cycles reached before convergence\n")
  }

  stats <- as.numeric(stats)
  stats <- matrix(stats, cycle, 2*k+1, byrow=TRUE)
  colnames(stats) <- c("llik",paste("s^2_", Vcoef.names, sep=""), paste("der_", Vcoef.names, sep=""))

  WX <- W %*% X
  XtWX <- crossprod(X,WX)
  cov <- XtWX
  cov <- solve(cov, cbind(t(WX),diag(1,dim(XtWX)[1])))
  beta.cov <- matrix(cov[,(dim(t(WX))[2]+1):dim(cov)[2]],dim(X)[2],dim(X)[2])
  cov <- matrix(cov[,1:dim(t(WX))[2]],dim(X)[2],dim(X)[1])

  beta <- cov %*% y
  beta <- matrix(beta,length(beta),1)
  row.names(beta) <- Xcolnames
  beta.se <- sqrt(abs(diag(beta.cov)))
  pos.cov <- (diag(beta.cov)<0)
  beta.se[pos.cov] <- NA
  beta.se <- matrix(beta.se,length(beta.se),1)
  row.names(beta.se) <- Xcolnames
  rms <- rss/rankQ

  fitted.values <- X%*%beta
  Q <- In -  X %*% cov

  predicted <- NULL
  predictedVariance <- NULL
  predictedVariance2 <- NULL
  if(identity) {
      gam <- sigma[k]  ## coefficient of identity, last variance term
      if(SWsolveINDICATOR) {
          ## Sigma is undefined
          Sigma <- 0
          for(i in 1:k)
          {
              Sigma <- Sigma + V[[i]] * sigma[i]
          }
      }
      predicted <- fitted.values + (Sigma - gam*In) %*% W%*%(y - fitted.values)
      ## predictedVariance <- Sigma - (Sigma - gam*In) %*% W %*% (Sigma - gam*In)
      ## really the last term should be transposed but they are square symmetric matrices here so....
      ## If you multiply it out and take the diagonal only it simplifies to....
      predictedVariance <- 2*gam - gam^2*diag(W)
      ## Variance of a new observation conditional on data (assuming
      ## known beta)
      predictedVariance2 <- diag(gam^2 * WX %*% beta.cov %*% t(WX))
      ## Additional variance of a new observation conditional on data taking into
      ## account variation in beta
  }

  ## scale dictated by pos
  sigma.cov <- (A.svd[1:k, 1:k] * 2)

  ## Last Step:  June 17th 2005
  ## Convert the estimates for the variance parameters, their standard
  ## errors etc to the usual scale

  FI <- A/2

  ## convert FI using pos
  FI.c <- matrix(0,dim(FI)[1],dim(FI)[2])
  FI.c <- FI / tcrossprod((sigma-1)*pos+1)
  ##for(i in 1:dim(FI)[1])
  ##  for(j in 1:dim(FI)[2])
  ##    FI.c[i,j] <- FI[i,j]/(((sigma[i]-1)*pos[i]+1)*((sigma[j]-1)*pos[j]+1))

  names(sigma) <- Vcoef.names
  sigma.cov <- try(ginv(FI.c),silent=TRUE)
  error1 <- (class(sigma.cov)=="try-error")
  if(error1) {
    cat("Warning: solution lies on the boundary; check sigma & pos\nNo standard errors for variance components returned\n")
    sigma.cov <- matrix(NA,k,k)
  }
  rownames(sigma.cov) <- colnames(sigma.cov) <- Vcoef.names

  ## Additional warning if any pos entries are TRUE and the corresponding term is close to zero
  ## Only give this warning it I haven't given the previous one.
  if(!error1) {
    if(any(sigma[pos]^2 < 1e-4)) cat("Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid\n")
  }
  
  result <- list(trace=stats, llik=llik, cycle=cycle, rdf=rankQ,
                 beta=beta, beta.cov=beta.cov, beta.se=beta.se,
                 sigma=sigma[1:k], sigma.cov=sigma.cov[1:k,1:k], W=W,
                 Q=Q, fitted=fitted.values, predicted=predicted,
                 predictedVariance=predictedVariance,
                 predictedVariance2=predictedVariance2,pos=pos,
                 Vnames=Vcoef.names, formula=formula,
                 Vformula=Vformula, Kcolnames=Kcolnames,
                 model=model,Z=Z, X=X, Sigma=Sigma)
  class(result) <- "regress"
  return(result)
}
