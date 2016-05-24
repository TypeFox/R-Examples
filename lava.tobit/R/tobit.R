###{{{ Objective function

##' @export
tobit_method.lvm <- "nlminb1"

##' @export
tobit_objective.lvm <- function(x,p,data,weight,indiv=FALSE,
                                algorithm=lava.options()$tobitAlgorithm,
                                seed=lava.options()$tobitseed,...) {
  if (!is.null(seed)) {
    if (!exists(".Random.seed")) runif(1)
    save.seed <- .Random.seed
    set.seed(seed)
  }
  require("mvtnorm")
  zz <- manifest(x) 
  d <- as.matrix(rbind(data)[,zz,drop=FALSE]);
  colnames(d) <- zz  
  yy <- endogenous(x)
  yy.w <- intersect(yy,colnames(weight))
  yy.idx <- match(yy.w,zz)
  Status <- matrix(0,ncol=length(zz),nrow=nrow(d))

  Status[,yy.idx]<- as.matrix(weight)[,yy.w,drop=FALSE]
  patterns <- unique(Status,MARGIN=1)
  cens.type <- apply(Status,1,
                     function(x) which(apply(patterns,1,function(y) identical(x,y))))  
  mp <- modelVar(x,p,data=as.data.frame(d)) 
  Sigma <- mp$C ## Model specific covariance matrix
  xi <- mp$xi ## Model specific mean-vector
  val <- c()
##  browser()
  for (i in 1:nrow(patterns)) {    
    ## Usual marginal likelihood for status==1
    pat <- patterns[i,]
    idx <- which(cens.type==i)
    noncens.idx <- which(pat==0)
    left.cens.idx <- which(pat==-1)
    right.cens.idx <- which(pat==1)
    cens.idx <- which(pat!=0)
    cens.which.left <- match(left.cens.idx,cens.idx)
    cens.which.right <- match(right.cens.idx,cens.idx)  
    noncens.y <- zz[noncens.idx]
    left.cens.y <- zz[left.cens.idx]
    right.cens.y <- zz[right.cens.idx] ##setdiff(zz,noncens.y)
    y <- d[idx,,drop=FALSE]
    val1 <- val0 <- 0;
    ## browser()
    if (length(noncens.y)>0) {
      ## p(y), using: int[p(y,y*)]dy* =  p(y) int[p(y*|y)]dy*
      val1 <- dmvnorm(d[idx,noncens.y,drop=FALSE], 
                      mean=xi[noncens.idx],
                      sigma=Sigma[noncens.idx,noncens.idx,drop=FALSE],
                      log=TRUE)
    }

    if (length(cens.idx)>0) {
      L <- diag(length(cens.idx))
      L[cbind(cens.which.left,cens.which.left)] <- (-1)
      ##  L[cens.which.left,cens.which.left] <- (-1)
      if (length(noncens.y)==0) {
        val0 <- sapply(idx,
                            function(ii)
                            log(pmvnorm(lower=as.numeric(L%*%d[ii,cens.idx]),
                                                     mean=as.numeric(L%*%xi),
                                        sigma=L%*%Sigma%*%L,algorithm=algorithm)))
      } else {
         M <- mom.cens(x,p,data=y,cens.idx,conditional=TRUE,deriv=FALSE)
         val0 <- c()
         for (j in 1:length(idx)) {           
           val0 <- c(val0,log(pmvnorm(lower=as.numeric(L%*%y[j,cens.idx]),mean=as.numeric(L%*%M$mu.censIobs[j,]),sigma=L%*%M$S.censIobs%*%L,algorithm=algorithm)) )
         }        
        }
    }    
    val <- c(val,-val1-val0)
  }
  if (!is.null(seed))
    .Random.seed <<- save.seed

  
  if (!indiv)
    return(sum(val))
  val
}

###}}} Objective

###{{{ Gradient & Hessian

##' @export
tobit_gradient.lvm <- function(x,p,data,weight,weight2=NULL,indiv=FALSE,
                               algorithm=lava.options()$tobitAlgorithm,
                               seed=lava.options()$tobitseed,...) {

  if (!is.null(seed)) {
    if (!exists(".Random.seed")) runif(1)
    save.seed <- .Random.seed
    set.seed(seed)
  }
  require("mvtnorm")
  zz <- manifest(x)
  d <- as.matrix(data[,zz,drop=FALSE]); colnames(d) <- zz
  yy <- endogenous(x)
  yy.w <- intersect(yy,colnames(weight))
  yy.idx <- match(yy.w,zz)
  Status <- matrix(0,ncol=length(zz),nrow=nrow(d))
  Status[,yy.idx]<- as.matrix(weight)[,yy.w,drop=FALSE]

  W0 <- weight2
  if (!is.null(W0)) {
    yy.w2 <- intersect(yy,colnames(weight2))
    yy.idx2 <- match(yy.w2,zz)    
    W0 <- matrix(1,ncol=length(zz),nrow=nrow(d))
    W0[,yy.idx2] <- weight2[,yy.w2,drop=FALSE]
    colnames(W0) <- zz
  }
##  yy <- endogenous(x)
##  yy.idx <- match(yy,zz)
##  Status <- matrix(0,ncol=length(zz),nrow=nrow(d))
##  Status[,yy.idx] <- as.matrix(weight)[,yy,drop=FALSE]
  patterns <- unique(Status,MARGIN=1)
  cens.type <- apply(Status,1,
                     function(x) which(apply(patterns,1,function(y) identical(x,y))))  
  val <- 0
##  score <- c()
  score <- matrix(ncol=length(p),nrow=nrow(data))
  for (i in 1:nrow(patterns)) {
    ## Usual marginal likelihood for status==1
    pat <- patterns[i,]
    idx <- which(cens.type==i)
    noncens.idx <- which(pat==0)
    left.cens.idx <- which(pat==-1)
    right.cens.idx <- which(pat==1)
    cens.idx <- which(pat!=0)
    cens.which.left <- match(left.cens.idx,cens.idx)
    cens.which.right <- match(right.cens.idx,cens.idx)  
    noncens.y <- zz[noncens.idx]
    left.cens.y <- zz[left.cens.idx]
    right.cens.y <- zz[right.cens.idx] ##setdiff(zz,noncens.y)
    y <- d[idx,,drop=FALSE]
    w <- W0[idx,,drop=FALSE]
    dummy <- cens.score(x,p,data=y,cens.idx=cens.idx,cens.which.left=cens.which.left, algorithm=algorithm,weight=w)
    score[idx,] <- dummy   
    ##    browser()
    ##    y.pat <- unique(y,MARGIN=1)
    ##    system.time(
    ##    y.type <- apply(y,1,
    ##                    function(x) which(apply(y.pat,1,function
    ##                                            (z) identical(x,z))))
    ##                )    
    ##    print(system.time(
    ##                      dummy <- cens.score(x,p,data=y.pat,cens.idx=cens.idx, cens.which.left=cens.which.left, algorithm=algorithm)
    ##                      ))
    ##    score <- rbind(score,dummy[y.type,,drop=FALSE]) 
    ##    score <- rbind(0)
    ##    browser()
    ##    dummy <- cens.score(x,p,data=y,cens.idx=cens.idx, cens.which.left=cens.which.left)
    ##    score <- rbind(score,dummy)
  }
  
  if (!is.null(seed))
    .Random.seed <<- save.seed
  if (indiv)
    return(-score)
  return(-colSums(score))
}

##' @export
tobit_hessian.lvm <- function(x,p,data,weight,...) {
  S <- -tobit_gradient.lvm(x,p=p,data=data,weight=weight,indiv=TRUE,...)
  J <- t(S)%*%S
  attributes(J)$grad <- colSums(S)
  return(J)  
}

###}}} Gradient & hessian

###{{{ Log-likelihood

##' @export
tobit_logLik.lvm <- function(object,p,data,weight,...) {
  res <- -tobit_objective.lvm(x=object,p=p,data=data,weight=weight,...)
  args <- list(...)
  args$object <- object; args$p <- p; args$data <- data; args$type <- "exo"; args$weight <- NULL
  ##xl <- lava:::gaussian_logLik.lvm(object,p=p,data=data,type="exo",weight=NULL,...)
  res <- res - do.call(gaussian_logLik.lvm,args)
  return(res)
}

###}}} Log-likelihood

###{{{ score for fixed censoring pattern

cens.score <- function(x,p,data,cens.idx,cens.which.left,weight,...) {
  obs.idx <- setdiff(1:NCOL(data),cens.idx)
  n <- NROW(data)
##  print(system.time(
  M <- mom.cens(x,p,data=data,cens.idx=cens.idx,conditional=TRUE,deriv=TRUE)
##))
  ## Censored part:
##  browser()
  if (length(cens.idx)>0) {

    combcens <- 1*(length(cens.which.left)>0) +
      -1*(length(cens.which.left)<length(cens.idx))
##    browser()
    ## 0: left and right, -1: left only, 1: right only
 
##     mu <- M$mu.cens
##     S <- M$S.cens
##     dmu <- M$dmu.cens
##     dS <- M$dS.cens
    z <- matrix(data[,cens.idx],nrow=n)
    ##    z <- data[,cens.idx,drop=FALSE]
    DCDFs <- c()

    S <- M$S.censIobs
    dS <- M$dS.censIobs

    if (combcens==0) {
      L <- -diag(length(cens.idx))
      L[cbind(cens.which.left,cens.which.left)] <- 1
      z <- z%*%L
      S <- L%*%S%*%L
      dS <- (L%x%L)%*%dS
    }
    w <- NULL
    for (i in 1:n) {
      mu <- M$mu.censIobs[i,]
      dmu <- matrix(M$dmu.censIobs[,,i],nrow=length(cens.idx))
      if (combcens==0) {
        mu <- L%*%mu
        dmu <- L%*%dmu
      }
      zi <- z[i,]
      if (!is.null(weight))
        w <- diag(weight[i,cens.idx],nrow=length(cens.idx))
      if (combcens==-1) {
        zi <- -zi
        mu <- -mu
        dmu <- -dmu
      }
      yy <- c(1,1); yy[cens.which.left] <- 0
##      print(yy)
##      print(mu)      
##      print(system.time(
      DCDF <- Dthetapmvnorm(zi,
                            ##        mu=L%*%M$mu.censIobs[i,],
                            mu=mu,
                            ##                            S=L%*%M$S.censIobs%*%L,
                            S=S,
                            ##                            dS=(L%x%L)%*%M$dS.censIobs,
                            dS=dS,
                            ##                            dmu=L%*%matrix(M$dmu.censIobs[,,i],nrow=length(cens.idx)),
                            dmu=dmu,
                            weight=w,
                            ...)
      alpha <- attributes(DCDF)$cdf
      DCDFs <- rbind(DCDFs, 1/alpha*DCDF)
    }
    ##    DCDF <- Dthetapmvnorm(z,mu=M$mu.cens,S=S,dmu=M$dmu.cens,dS=M$dS.cens)
    ##    S0 <- 1/attributes(DCDF)$cdf*DCDF
    S0 <- DCDFs    
  } else S0 <- 0

  
  ## Observed part:
  ##  S1 <- matrix() ...
  if (length(obs.idx)>0) {
    y1 <- as.matrix(data[,obs.idx,drop=FALSE])
    if (n==1) y. <- y1 else y. <- colMeans(y1)
    mu <- M$mu.obs
    S <- M$S.obs
    iS <- Inverse(S)
    dS <- M$dS.obs
    dmu <- M$dmu.obs
    S1 <- c()
    part0 <- -1/2*as.vector(iS)%*%dS    
    for (i in 1:n) {
      z <- as.numeric(y1[i,])
      u <- z-mu
##      browser()
      if (!is.null(weight)) {
        W <- diag(weight[i,obs.idx],nrow=length(obs.idx))
        S1 <- rbind(S1,
                    as.numeric(crossprod(W%*%u,iS)%*%dmu -
                               1/2*as.vector((iS-iS%*%tcrossprod(u)%*%iS)%*%W)%*%dS))
      } else {      
        S1 <- rbind(S1,
                    as.numeric(part0 + crossprod(u,iS)%*%dmu +
                               1/2*as.vector(iS%*%tcrossprod(u)%*%iS)%*%dS))
      }
    }
  } else S1 <- 0
##   cat(rep("*",50),sep="")
##   print(S1)
##   cat(rep("-",50),sep="")
##   print(S0)
##   cat(rep("#",50),sep="")
##   print(dim(data))
  return((S1+S0))
}

###}}} score for fixed censoring pattern

###{{{ Derivatives of normal CDF

## Calculates first and second order partial derivatives of normal CDF
Dpmvnorm <- function(Y,S,mu=rep(0,NROW(S)),std=FALSE,seed=lava.options()$tobitseed,
                     algorithm=lava.options()$tobitAlgorithm,
                     ...) {
##  browser()
  k <- NROW(S)
  if (!is.null(seed) & k>1) {
    if (!exists(".Random.seed")) runif(1)
    save.seed <- .Random.seed
  }
  require("mvtnorm")
  if (!std) {
    L <- diag(S)^0.5
    Li <- diag(1/L,NROW(S))
    L <- diag(L,NROW(S))
    S0 <- S
    S <- Li%*%S%*%Li
    ##    Y0 <- apply(Y,1,function(x)-mu)
    Y <- (Li%*%(Y-mu))
  }
  ##  Y <- as.numeric(Y)
  Y <- as.vector(Y)
  
  if (k==1) {
    D <- dnorm(Y,sd=as.vector(S)^0.5)
    H <- -Y*D
    return(list(grad=D,hessian=H))
  }
  ## Cond. var of Y[-j] given Y[j]
  D <- numeric(k)
  for (j in 1:k) {
##    tcrossprod(S[-j,j,drop=FALSE])
    Sj <- S[-j,-j,drop=FALSE] - tcrossprod(S[-j,j,drop=FALSE]) ##/S[j,j] S=correlation
    muj <- Y[-j] - S[-j,j,drop=FALSE]*Y[j]
#    set.seed(seed)
    D[j] <- dnorm(Y[j])*pmvnorm(upper=as.vector(muj),sigma=Sj,algorithm=algorithm)
  }
  H <- matrix(0,k,k) 
  if (k<3) {
    H[1,2] <- H[2,1] <- dmvnorm(unlist(Y),sigma=S)
    diag(H) <- -Y*D -H[1,2]*S[1,2]
    ##as.vector((H*S)%*%rep(1,k))
  } else {
    ##    H[] <- 0
    phis <- Phis <- H
    for (i in 1:(k-1)) {
      for (j in (i+1):k) {
        Snij <- S[-c(i,j),c(i,j),drop=FALSE]
        B <- Snij%*%Inverse(S[c(i,j),c(i,j)])
        Sij <- S[-c(i,j),-c(i,j),drop=FALSE] - B%*%t(Snij)
        muij <- Y[-c(i,j)] - B%*%Y[c(i,j)]
#        set.seed(seed)
        Phis[i,j] <- Phis[j,i] <- pmvnorm(upper=as.vector(muij),sigma=Sij,algorithm=algorithm)
        phis[i,j] <- phis[j,i] <- dmvnorm(Y[c(i,j)],sigma=S[c(i,j),c(i,j)])
      }
    }
    H <- Phis*phis
    diag(H) <- -Y*D - as.vector((S*H)%*%rep(1,k))    
  }
  if (!std) {
    if (!is.null(seed))
      .Random.seed <<- save.seed
    a <- pmvnorm(upper=Y,mean=as.numeric(mu),sigma=S,algorithm=algorithm)
    return(list(grad=Li%*%D, hessian=Li%*%H%*%Li, R=S, CDF=a, S=S0, mu=mu, L=L, Li=Li))
  }
  if (!is.null(seed))
    .Random.seed <<- save.seed
  return(list(grad=D,hessian=H))
}

## Calculates first and second order partial derivatives of normal CDF
## w.r.t. parameter-vector!
Dthetapmvnorm <- function(yy,mu,S,dmu,dS,seed=lava.options()$tobitseed,
                          algorithm=lava.options()$tobitAlgorithm, weight,
                          ...) {
  if (!is.null(seed)) {
    if (!exists(".Random.seed")) runif(1)
    save.seed <- .Random.seed
  }
  ##yy <- as.matrix(yy)
  ##  pp <- modelPar(x,p)
  ##  M <- moments(x,p)
  ##  mu <- M$xi
  ##  S <- M$C
  iS <- Inverse(S)  
  ##  DCDF <- Dpmvnorm(y,S,mu)
  L <- diag(S)^0.5
  Li <- diag(1/L,NROW(S))
  L <- diag(L,NROW(S))
  R <- Li%*%S%*%Li
  LR <- L%*%R ## = S%*%Li
  if (!is.null(weight)) {    
    K1 <- -0.5*t(dS)%*%cbind(as.vector(iS%*%weight))
    K2 <- 0.5*t(dS)%*%(iS%x%iS)
    K3 <- t(dmu)%*%(iS%*%weight)
  } else {
    K1 <- -0.5*t(dS)%*%cbind(as.vector(iS))
    K2 <- 0.5*t(dS)%*%(iS%x%iS)
    K3 <- t(dmu)%*%(iS)
  }
  S0 <- function(y) {
    z <- Li%*%(y-mu)
#    set.seed(seed)
    a <- pmvnorm(upper=y,mean=as.numeric(mu),sigma=S,algorithm=algorithm)
    DC <- Dpmvnorm(z,R,std=TRUE,algorithm=algorithm)
    MM <- -LR%*%(DC$grad)
    VV <- LR%*%(DC$hessian)%*%t(LR) + a*S
##    message("DC")
##    print(DC$hessian)
    part1 <- K1*a
    if (!is.null(weight)) {
      VV <- VV%*%weight
    }
    part2 <- K2%*%as.vector(VV)
    part3 <- K3%*%as.vector(MM)
    res <- part1+part2+part3
    return(c(a,res))
  }
##  DD <- t(apply(yy,1,S0))
  DD <- matrix(S0(as.numeric(yy)),nrow=1)
  a <- DD[,1]
  res <- DD[,-1,drop=FALSE]
  attributes(res)$cdf <- a
  if (!is.null(seed))
    .Random.seed <<- save.seed
  return(res)
}

###}}} Derivatives of normal CDF

###{{{ Derivatives of Marginal and Conditional moments of a normal distribution 

mom.cens <- function(x,p,cens.idx,data,deriv=TRUE,conditional=TRUE,right=TRUE,...) {
  obs.idx <- setdiff(1:NCOL(data),cens.idx)
##  browser()
  M <- moments(x,p,data=as.data.frame(data))
  if (deriv)
    D <- deriv(x,p=p,mom=M,meanpar=TRUE) ##,mu=colMeans(data))

  if (length(cens.idx)<1) {
    res <- list(S.obs=M$C, mu.obs=M$xi, S.cens=NULL, mu.cens=NULL,
                S.censIobs=NULL, mu.censIobs=NULL)
    if (deriv)
      res <- c(res, list(dS.obs=D$dS, dmu.obs=D$dxi, dS.cens=NULL, dmu.cens=NULL,
                         dS.censIobs=NULL, dmu.censIobs=NULL
                         ))
    return(res)
  }
  if (length(obs.idx)<1) {
    res <- list(S.obs=NULL, mu.obs=NULL, S.cens=M$C, mu.cens=M$xi,
                S.censIobs=M$C, mu.censIobs=matrix(M$xi,ncol=length(M$xi),nrow=NROW(data),byrow=TRUE)
                )
    if (deriv)
      res <- c(res, list(dS.obs=NULL, dmu.obs=NULL, dS.cens=D$dS, dmu.cens=D$dxi,
                         dS.censIobs=D$dS, dmu.censIobs=array(D$dxi,dim=c(dim(D$dxi),NROW(data)))
               ))
    return(res)
  }
  S.cens <- M$C[cens.idx,cens.idx,drop=FALSE]
  S.obs <- M$C[obs.idx,obs.idx,drop=FALSE]
  S.censobs <- M$C[cens.idx,obs.idx,drop=FALSE]
  S.obscens <- M$C[obs.idx,cens.idx,drop=FALSE]
  mu.obs <- M$xi[obs.idx]
  mu.cens <- M$xi[cens.idx]
  iS.obs <- Inverse(S.obs)

  S01iS0 <- S.censobs%*%iS.obs
  S.cond <- S.cens - S01iS0%*%S.obscens
  y0.obs <- apply(data,1,function(z) z[obs.idx]-mu.obs)
  Zy0 <- S01iS0%*%y0.obs
  mu.cond <- t(rbind(apply(Zy0,2,function(z) mu.cens+z)))
  
  res <- list(S.obs=S.obs, mu.obs=mu.obs, S.cens=S.cens, mu.cens=mu.cens)
  if (conditional)
    res <- c(res, list(S.censIobs=S.cond, mu.censIobs=mu.cond))
  if (!deriv)
    return(res)
           
  M.idx <- matrix(1:nrow(M$C)^2,nrow(M$C))
  Cens <- as.vector(M.idx[cens.idx,cens.idx])
  Obs <- as.vector(M.idx[obs.idx,obs.idx])
  CensObs <- as.vector(M.idx[cens.idx,obs.idx])
  ObsCens <- as.vector(M.idx[obs.idx,cens.idx])
  dS.cens <- D$dS[Cens,,drop=FALSE]
  dS.obs <- D$dS[Obs,,drop=FALSE]
  dS.censobs <- D$dS[CensObs,,drop=FALSE]
  dS.obscens <- D$dS[ObsCens,,drop=FALSE]

  res <- c(res, list(##dT.obs=D$dT[Obs,], 
                     dS.obs=dS.obs, dmu.obs=D$dxi[obs.idx,,drop=FALSE],
                     dS.cens=dS.cens, dmu.cens=D$dxi[cens.idx,,drop=FALSE])
           )
  if(!conditional)
    return(res)

  K0 <- S.censobs%*%iS.obs
  K1 <- K0%x%diag(length(cens.idx))
  K2 <- K0%x%K0
###  K3 <- diag(length(obs.idx))%x%K0
  K3 <- diag(length(cens.idx))%x%K0
  dSc <-
    dS.cens -
      K1%*%dS.censobs +
        K2%*%dS.obs -
          K3%*%dS.obscens
  dxic.1 <- D$dxi[cens.idx,] - K0%*%D$dxi[obs.idx,]

  dxic <- array(dxic.1, dim=c(nrow(dxic.1),ncol(dxic.1),nrow(data)))
  for (k in 1:length(p)) {
    dxic[,k,] <- dxic[,k,] +
      (
       matrix(dS.censobs[,k],nrow=length(cens.idx),ncol=length(obs.idx)) -
       K0%*%matrix(dS.obs[,k],ncol=length(obs.idx))
       )%*%iS.obs%*%y0.obs
  }

  ##  y <- as.matrix(d[1,obs.idx])
  ##  y0 <- y-M$xi[obs.idx]  
  ##  Ky <- (y0)%*%solve(S.obs)
  ##  K1 <- Ky%x%diag(length(cens.idx))
  ##  K2 <- Ky%x%K0
  ##  dxic.2 <- K1%*%dS.censobs - K2%*%dS.obs
  ##  dxic <- dxic.1+dxic.2
  
  res <- c(res, list(#D=invisible(D),
                     dS.censIobs=dSc, dmu.censIobs=dxic))
  return(res)
}

###}}} Derivatives of Marginal and Conditional moments of a normal distribution 
