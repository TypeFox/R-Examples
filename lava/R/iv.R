###{{{ Objective

IV_method.lvm <- NULL
IV_objective.lvm <- function(x,p,data,...) {
  IV2(x,data,...)
}
IV_variance.lvm <- function(x,p,data,opt,...) {
  opt$vcov
}

IV0_objective.lvm <- function(x,p,data,...) {
  IV(x,data)
}
IV0_variance.lvm <- function(x,p,data,opt,...) {
  opt$vcov
}

###}}} Objective

CondVar <- function(S,idx) {
  idx2 <- setdiff(seq_len(ncol(S)),idx)
  S11 <- S[idx2,idx2];
  S22 <- S[idx,idx]
  S12 <- S[idx2,idx]
  S11-S12%*%solve(S22)%*%t(S12)
}

varest <- function(x,data) {
  p <- IV(x,data)$estimate
  idx <- match(names(p),coef(x,mean=TRUE))
  x0 <- parfix(Model(x),idx,p)
  index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)

  A <- t(index(x)$A)
  Afix <- A; Afix[t(index(x)$M0)==1] <- 0
  A[A!=0] <- 1
  k <- nrow(A)
  I <- diag(nrow=k)
  Ap <- modelVar(x)$A ## Estimated parameter matrix

  indicators <- setdiff(vars(x)[rowSums(A)==1],exogenous(x))
  responses <- endogenous(x,top=TRUE)
  y.indicators <- responses[rowSums(A[responses,])==1]
  Sigma <- var(data[,manifest(x)])

  var.eta <- c()
  for (eta in latent(x)) {
    m.sub <- subset(Model(x),c(eta,indicators))
    reachable <- acc(x$M,eta)
    ys <- intersect(names(reachable),y.indicators)
    lambdas <- c()
    for (y in ys) {
      pp <- path(Model(x), from=eta, to=y)
      lambda1 <- 0
      for (i in seq_along(pp)) {
        lambda <- 1
        for (j in seq_len(length(pp[[i]])-1))
          lambda <- lambda*Ap[pp[[i]][j],pp[[i]][j+1]]
        lambda1 <- lambda1+lambda
      }
      lambdas <- c(lambdas,lambda1)
    }
    val <- outer(1/lambdas,1/lambdas)*Sigma[ys,ys]
    var.eta <- c(var.eta, mean(val[upper.tri(val)]))
  }

  S <- rep(0,k); S[match(manifest(x),vars(x))] <- diag(Sigma); S[match(latent(x),vars(x))] <- var.eta; names(S) <- vars(x)
  I <- diag(nrow=k)
  IA <- (I-t(Ap))
  IA%*%cbind(S)%*%t(IA)

}


## Instrumental Variable Estimator / 2SLS
##' @export
IV <- function(m,data,R2thres=0,...) {
  if (length(constrain(m))>0) stop("Nonlinear constrains not supported!")
  if (inherits(m,"lvmfit")) {
      m <- Model(m)
  }
  R2 <- cor(data[,manifest(m)])^2

  A <- t(index(m)$A)
  Afix <- A; Afix[t(index(m)$M0)==1] <- 0
  A[A!=0] <- 1
  P <- index(m)$P
  k <- nrow(A)
  I <- diag(nrow=k)
  B <- rbind(I,solve(I-A))
  VV <- B%*%P%*%t(B)
  u.var <- index(m)$vars
  all.idx <- seq_along(u.var)
  lat.idx <- with(index(m), which(vars%in%latent))
  if (length(lat.idx)==0) stop("Estimator only defined for models with latent variable")
  y.var <- endogenous(m)
  y.idx <- which(index(m)$vars%in%y.var)
  x.idx <- which(vars(m)%in%exogenous(m))

  ## Set of Indicator variables:
  indicators <- c()
  for (i in seq_len(nrow(A))) {
    ifix <- (Afix[i,]==1)
    if ((sum(ifix)==1) &  all(A[i,-which(ifix)]==0))
      indicators <- c(indicators, i)
  }
  y.indicators <- intersect(indicators,y.idx)

  y.scale <- list()
  for (eta in lat.idx) {
    pred.eta <- intersect(y.idx, which(Afix[,eta]==1)) ## Candidates for
    ## eta = y-epsilon
    if (length(pred.eta)<1)
      pred.eta <- intersect(lat.idx, which(Afix[,eta]==1))
    myidx <- c()
    for (y in pred.eta) {
      y.pred <- setdiff(eta,which(A[y,]==1)) ## No other variables predicting y?
      if (length(y.pred)==0)
        myidx <- c(myidx,y)
    }
    y.scale <- c(y.scale, list(myidx))
  }

  if (any(unlist(lapply(y.scale, function(x) length(x)))<1)) stop("At least one scale-measurement pr. latent variable")

  vv <- setdiff(seq_len(k),c(unlist(y.scale),x.idx))

  Ecor <- list()
  eta.surrogate <- c()
  latno <- 0
  for (e in lat.idx) {
    latno <- latno+1
    y0 <- y.scale[[latno]][1]
    if (!(y0%in%lat.idx)) {
      eta.surrogate <- c(eta.surrogate,vars(m)[y0])
      Ecor <- c(Ecor,list(y0))
    }
    else {
      m.sub <- subset(m,vars(m)[c(e,indicators)])
      i <- 0
      while (i<length(y.indicators)) {
        i <- i+1
        pp <- path(m.sub,from=vars(m)[e],to=vars(m)[y.indicators[i]])[[1]]
        if (!is.null(pp)) {
          Ecor <- c(Ecor,
                    list(which(vars(m)%in%pp[-1])))
          eta.surrogate <- c(eta.surrogate, tail(pp,1))
        }
      }
    }
  };
  names(eta.surrogate) <- latent(m)

  dd <- list()
  ll  <- list()
  coefname <- coef(m,mean=TRUE)
  mycoef <- rep(0,length(coefname))
  A0 <- A
  P0 <- P
  D <- c()
  V. <- list()
  Z. <- list()
  Y. <- list()
  count <- 0
  ff <- list()
  instruments <- c()
  parname <- c()
  for (v in vv) {
    pred <- which(A[v,]==1)
    if (sum(pred)>0) {
      Debug(vars(m)[v])
      pred.lat <- intersect(pred,lat.idx) # Any latent predictors?
      lpos <- match(v,lat.idx)
      lppos <- match(pred.lat,lat.idx)
      ecor <- c(v,unlist(Ecor[lppos]))
      if (!is.na(lpos)) {
        v0 <- match(eta.surrogate[lpos],vars(m))
        ecor <- c(ecor,Ecor[[lpos]])
      } else {
        v0 <- v
      }

      ecor <- unique(c(v0,ecor))
      XX <- vars(m)[A[v,]==1]
      intpred <- exogenous(m)
      newf <- c()
      if (length(pred.lat)>0) {
        intpred <- vars(m)
        for (i in seq_along(pred.lat)) {
          uncor <- which(colSums(VV[ecor,k+seq_len(k),drop=FALSE])==0)
          uncor <- setdiff(uncor,c(lat.idx))
          mypred <- vars(m)[uncor]
          XX[XX==vars(m)[pred.lat[i]]] <- eta.surrogate[lppos[i]]
          ##          allpred <- c(allpred, mypred)
          intpred <- intersect(intpred,mypred)
          f <- toformula(eta.surrogate[lppos[i]],mypred)
          ff <- c(ff,
                  f)
          f2 <- list(f)
          names(f2) <- vars(m)[i]
          newf <- c(newf,f2)
        }
      }

      intpred <- intersect(intpred,manifest(m))
      R2max <- apply(R2[XX,intpred,drop=FALSE],2,max)
      if (any(R2max<R2thres)) intpred <- intpred[R2max>=R2thres]
      newf <- list(intpred); names(newf) <- vars(m)[v]
      instruments <- c(instruments, newf)
      covariates <- unique(c(setdiff(colnames(A)[A[v,]==1],latent(m)),intpred))##allpred)
      if (length(covariates)==0) stop("No instruments")
      Z <- model.matrix(toformula("",c("1",XX)),data)
      Y <- as.matrix(data[,vars(m)[v0]])
      V <- model.matrix(toformula("",c("1",unique(covariates))),data)
      count <- count+1
      V. <- c(V.,list(V))
      Z. <- c(Z.,list(Z))
      Y. <- c(Y.,list(Y))
      XX <- vars(m)[A[v,]==1 & Afix[v,]!=1]
      parname <- c(parname, c(vars(m)[v0],paste(vars(m)[v],XX,sep=lava.options()$symbol[1])))
    } else {
      if (vars(m)[v]%in%latent(m)) {
        lpos <- match(v,lat.idx)
        v0 <- match(eta.surrogate[lpos],vars(m))
        Y <- matrix(data[,vars(m)[v0]],ncol=1)
        Y. <- c(Y.,list(Y))
        V. <- c(V.,list(cbind(rep(1,nrow(Y)))))
        Z. <- c(Z.,list(cbind(rep(1,nrow(Y)))))
        parname <- c(parname, names(eta.surrogate)[lpos])
       }
    }
  }

  LS <- function(X) {
    with(svd(X), v%*%diag(1/d,nrow=length(d))%*%t(u))
  }
  projection <- function(X) X%*%LS(X)
  P0 <- lapply(V.,LS)
  Zhat <- list(); for (i in seq_along(Z.)) Zhat <- c(Zhat, list(V.[[i]]%*%(P0[[i]]%*%Z.[[i]])))
  ZhatLS <- lapply(Zhat,LS)
  theta <- list(); for (i in seq_along(Y.)) theta <- c(theta, list(ZhatLS[[i]]%*%Y.[[i]]))
  u <- c()
  for (i in seq_along(Y.))
    u <- cbind(u, Y.[[i]]-Z.[[i]]%*%theta[[i]])
  covu <- crossprod(u)/nrow(u)

  theta.npar <- unlist(lapply(theta,length))
  theta.ncum <- c(0,cumsum(theta.npar))
  vartheta <- matrix(0,ncol=sum(theta.npar),nrow=sum(theta.npar))
  for (i in seq_along(theta)) {
    for (j in seq(i,length(theta))) {
      idx1 <- seq_len(theta.npar[i]) + theta.ncum[i]
      idx2 <- seq_len(theta.npar[j]) + theta.ncum[j]
      uZZ <- covu[i,j]* (ZhatLS[[i]]%*%t(ZhatLS[[j]]))
      vartheta[idx1,idx2] <- uZZ
      if (i!=j) {
        vartheta[idx2,idx1] <- t(uZZ)
      }
    }
  }


  parname[which(parname%in%eta.surrogate)] <- names(eta.surrogate)[which(eta.surrogate%in%parname)]

  coef <- cbind(unlist(theta),diag(vartheta)^0.5); rownames(coef) <- parname; colnames(coef) <- c("Estimate","Std.Err")
  res <- list(estimate=coef[,1], vcov=vartheta)
  attributes(res)$surrogates <- eta.surrogate
  attributes(res)$instruments <- instruments
  return(res)
}

IV2 <- function(m,data,control=list(),...) {
  if (is.null(control$R2thres)) control$R2thres <- 0
  res <- IV(m,data,R2thres=control$R2thres)
  p <- res$estimate
  idx <- match(names(p),coef(m,mean=TRUE))
  x0 <- parfix(m,idx,p)
  index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)
  idx0 <- order(idx)
  p0 <- p[idx0]
  V0 <- res$vcov[idx0,idx0]
  if (is.null(control$variance) || control$variance) {
    suppressWarnings(e0 <- estimate(x0,data,...,silent=TRUE,quick=TRUE))
    p0 <- c(p0,e0)
    V0 <- V0%++%matrix(0,ncol=length(e0),nrow=length(e0))
  }
  R2 <- noquote(formatC(cor(data[,manifest(m)])^2))
  colnames(R2) <- rownames(R2) <- manifest(m)
  l1 <- noquote(rbind(paste(latent(m),collapse=","),
                      paste(attributes(res)$surrogates,collapse=","),
                      ""))
  rownames(l1) <- c("Latent variables","Surrogate variables:","")
  colnames(l1) <- ""
  ii <- attributes(res)$instruments
  I <- noquote(matrix(NA,ncol=2,nrow=length(ii)))
  rownames(I) <- rep("",nrow(I))
  colnames(I) <- c("Response","Instruments")
  for (i in seq_along(ii)) {
    I[i,] <- c(names(ii)[i],paste(ii[[i]],collapse=","))
  }
  mymsg <- list(l1,I);
  list(estimate=p0,vcov=V0,summary.message=function(...)  {
       mymsg })
}
