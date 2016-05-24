###{{{ starter.multigroup

##' @export
starter.multigroup <- function(x, starterfun=startvalues2, meanstructure=TRUE,silent=TRUE,...) {
  ## Initial values:
  W <- c() ## Weight-vector
  s <- list()
  for (i in seq_len(x$ngroup)) {
    mydata <- x$data[[i]][,manifest(x$lvm[[i]]),drop=FALSE]
    W <- c(W, nrow(mydata))
    if (nrow(mydata)<3) {
      ii <- index(x$lvm[[i]])
      nn <- ifelse(meanstructure, ii$npar+ii$npar.mean, ii$npar)
      s0 <- rep(1,nn)
    }
    else {
      S <- x$samplestat[[i]]$S
      mu <- if (meanstructure) x$samplestat[[i]]$mu else NULL;
      ##      S <- cov(mydata); mu <- if (meanstructure) colMeans(mydata) else NULL;
      s0 <- starterfun(x$lvm[[i]], S=S, mu=mu,silent=TRUE)
    }
    s <- c(s, list(s0))
  }

  Wtotal <- sum(W); W <- W/Wtotal

  pg <- vector("list", x$npar); for (i in seq_len(length(pg))) pg[[i]] <- rep(0,x$ngroup)
  meang <- vector("list", x$npar.mean); for (i in seq_len(length(meang))) meang[[i]] <- rep(0,x$ngroup)

  for (i in seq_len(x$ngroup)) {
    pp <- modelPar(x$lvm[[i]],s[[i]])
    pos <- sapply(x$parlist[[i]], function(y) as.numeric(substr(y,2,nchar(y))))
    for (j in seq_len(length(pos)))
      pg[[ pos[j] ]][i] <-  pp$p[j]

    pos <- sapply(x$meanlist[[i]], function(y) as.numeric(substr(y,2,nchar(y))))
    ptype <- sapply(x$meanlist[[i]], function(y) substr(y,1,1)=="m")
    if (!(any(ptype)))
      pos <- NULL
    else
      pos <- pos[ptype]
    if (length(pos)>0)
    for (j in seq_len(length(pos))) {
      meang[[ pos[j] ]][i] <-  pp$meanpar[j]
    }
  }

  ## Weighted average
  wp <- unlist(lapply(pg, function(y) {
    ppos <- !is.na(y)
    myweight <- W[ppos]/sum(W[ppos])
    sum(y[ppos]*myweight)
  }))
  wmean <- unlist(lapply(meang, function(y) {
    ppos <- !is.na(y)
    myweight <- W[ppos]/sum(W[ppos])
    sum(y[ppos]*myweight)
  }))
  res <- c(wmean,wp)
  res[!is.finite(res) | is.nan(res) | is.na(res) | is.complex(res)] <- .5
  return(as.numeric(res))
}

###}}}

###{{{ startmean

startmean <- function(x,p,mu) {
  if (is.null(mu))
    return(p)
  meanpar <- numeric(index(x)$npar.mean)
  mymeans <- vars(x)[index(x)$v1==1]
  midx <- na.omit(match(names(mu),mymeans))
  meanpar[midx] <- mu[midx]
  AP <- matrices(x,p,meanpar)
  nu <- numeric(length(vars(x)))
  nu[vars(x)%in%manifest(x)] <- mu
  meanstart <- ((diag(nrow=nrow(AP$A))-t(AP$A))%*%nu)[index(x)$v1==1]
  names(meanstart) <- vars(x)[index(x)$v1==1]
  return( c(meanstart, p) )
}

###}}}

###{{{ startvalues3

`startvalues3` <-
function(x, S, debug=FALSE, tol=1e-6,...) {
  S <- reorderdata.lvm(x,S)
  if (nrow(S)!=length(manifest(x))) stop("Number of observed variables in data and models does not agree.")
  J <- index(x)$J ## Manifest selection
  P0 <- index(x)$P0 ## covariance 'adjacency'
  A <- t(index(x)$M) ## Adjacency matrix
  n <- nrow(S) ## Number of manifest variables
  m <- nrow(A) ## Number of variables
  A1 <- t(index(x)$M1) ## Adjacency matrix (without fixed parameters and duplicates)
  A0 <- t(index(x)$M0) ## Adjacency matrix (without fixed parameters)

  obs.idx <- index(x)$obs.idx;
  ##obs.idx <- as.vector(J%*%(seq_len(m)));
  latent.idx <- setdiff(seq_len(m), obs.idx)
  lat <- colnames(A)[latent.idx]

  exo.idx <- index(x)$exo.idx ## match(exogenous(x),vars(x))
  exo.idxObs <- index(x)$exo.obsidx ##match(exogenous(x),manifest(x))

  AP0 <- moments(x, rep(0,index(x)$npar))
  newP <- t(AP0$P)
  newA <- t(AP0$A)
  fixed <- t(x$fix)

  for (i in latent.idx) {
    fix.idx <- colnames(fixed)[which(!is.na(t(fixed[,i])))[1]]
    lambda0 <- newA[fix.idx,i]
    rel.idx <- which(A0[,i]==1)
    rel.all <- which(A[,i]==1)
    rel.pos <-  colnames(A)[rel.all]
    ## Estimation of lambda (latent -> endogenous)
    for (j in rel.idx) {
      lambda <- lambda0*S[fix.idx, j]/S[fix.idx,fix.idx]
      newA[j,i] <- lambda
    }
    lambdas <- newA[rel.pos,i]


    ## Estimation of beta (covariate -> latent)
    exo2latent <- which(A0[i,exo.idx]==1)
    exo.pos <- colnames(S)[exo.idxObs[exo2latent]]
    varX.eta <- S[exo.pos, exo.pos]
    InvvarX.eta <- Inverse(varX.eta,tol=1e-3)
    rel.pos <- setdiff(rel.pos,lat)
    covXY <- S[exo.pos, rel.pos,drop=FALSE]
    beta <- 0
    for (j in seq_len(length(rel.pos)))
      beta <- beta + 1/lambdas[j]*InvvarX.eta %*% covXY[,j]
    beta <- beta/length(rel.pos)

    for (k in seq_len(length(exo.pos))) {
      if (A0[i,exo.pos[k]]==1) {
        newA[i,exo.pos[k]] <- beta[k]
      }
    }

    beta.eta <- matrix(newA[i,exo.pos], ncol=1)

    ## Estimation of  zeta^2 (variance of latent variable)
    betavar <- matrix(beta.eta,nrow=1)%*%varX.eta%*%beta.eta

    zetas <- c()
    for (r1 in seq_len(length(rel.pos)-1))
      for (r2 in seq(r1+1,length(rel.pos))) {
        zetas <- c(zetas, S[rel.pos[r1], rel.pos[r2]]/ (lambdas[r1]*lambdas[r2]) - betavar)
      }
    zeta <- mean(zetas)

    newP[i,i] <- zeta
    for (j in rel.all) {
      pos <- colnames(newA)[j]
      vary <- S[pos,pos] - newA[pos,i]^2*(zeta+betavar)
      newP[pos,pos] <- ifelse(vary<0.25,0.25,vary)
    }

  }
  Debug(list("start=",start), debug)
  start <- pars(x, A=t(newA), P=newP)
  return(start)
}

###}}} startvalues3

###{{{ startvalues2


## Estimate sub-models (measurement models)
##' @export
`startvalues2` <-
  function(x, S, mu=NULL, debug=FALSE, silent=FALSE,...) {
    if (!silent) cat("Obtaining start values...\n")
    S <- reorderdata.lvm(x,S)
    ss <- startvalues(x,S)
    Debug(list("ss=",ss),debug);
    g <- measurement(x,silent=TRUE)
    keep <- c()
    if (length(g)>1) {
      for (i in seq_len(length(g))) {
        if (length(endogenous(g[[i]]))>2)
          keep <- c(keep,i)
      }
      g <- g[keep]
    }
    if (length(g)<2)
      return(startmean(x,ss,mu=mu))
    ## if (!silent) cat("Fitting marginal measurement models...\n")
    op <- options(warn=-1)
    e <- lapply(g, function(y) {
        estimate(y, data=list(S=S[manifest(y),manifest(y),drop=FALSE], mu=mu[manifest(y)], n=100), control=list(meanstructure=FALSE, starterfun="startvalues", estimator="Simple", method="nlminb1"), optcontrol=list(), debug=FALSE, silent=TRUE)
    })
    for (l in e) {
      ##    a <- coef(l$estimate)[,1]
      a <- coef(l)
      for (i in seq_len(length(a))) {
        pos <- match(names(a)[i],names(ss))
        if (!is.na(pos))
          ss[pos] <- a[i]
    }
    }
    options(op)
    startmean(x,ss,mu=mu)
  }

###}}} startvalues2

###{{{ startvalues0

##' @export
startvalues1 <- function(x,S,mu=NULL,tol=1e-6,delta=1e-6,...) {
    p0 <- startvalues(x,S,mu,...)
    p0[index(x)$npar.mean+variances(x)] <- 0.1
    p0[index(x)$npar.mean+offdiags(x)] <- 0
    p0
}

startvalues00 <- function(x,S,mu=NULL,tol=1e-6,delta=1e-6,...) {
    p0 <- startvalues(x,S,mu,...)
    p0 <- numeric(length(p0))
    P0 <- x$cov*1
    ##P0[!is.na(x$covfix)] <-
    ##P0 <- x$covfix; P0[is.na(P0)] <- 0
    ##diag(P0)[index(x)$endo.idx] <- diag(S)[index(x)$endo.obsidx]/2
    lu <- min(diag(P0)[index(x)$endo.idx])/2
    ##    diag(P0)[] <- 0.1
    ## diag(P0)[index(x)$endo.idx] <- 1
    diag(P0)[index(x)$eta.idx] <- 0.1 ##mean(diag(S)[index(x)$endo.idx])/2
    ee <- eigen(P0)
    tol <- 1e-6
    ii <- ee$values
    ii[ee$values<tol] <- tol
    P0 <- ee$vectors%*%diag(ii,nrow=length(ii))%*%t(ee$vectors)
    A0 <- P0
    pp <- pars(x,A=t(A0),P=P0,v=rep(0,index(x)$npar.mean))
    idx <- index(x)$npar.mean + c(offdiags(x),variances(x))
    p0[idx] <- pp[idx]
    return(p0)
}


##' @export
startvalues0 <- function(x,S,mu=NULL,tol=1e-6,delta=1e-6,...) {
    p0 <- startvalues(x,S,mu,...)
    A <- t(index(x)$M) ## Adjacency matrix
    P0 <- A0 <- matrix(0,nrow(A),ncol(A))
    A0[,index(x)$eta.idx] <- A[,index(x)$eta.idx]
    diag(P0)[index(x)$endo.idx] <- diag(S)[index(x)$endo.obsidx]/3
    ##lu <- min(diag(P0)[index(x)$endo.idx])
    lu <- 0.9
    diag(P0)[index(x)$eta.idx] <- lu##mean(diag(S)[index(x)$endo.idx])/2
    pp <- pars(x,A=t(A0),P=P0,v=rep(1,length(index(x)$vars)))
    nu <- numeric(length(vars(x)))
    pp[pp==1] <- p0[pp==1]
    ## if (!is.null(mu)) {
    ##     nu[vars(x)%in%manifest(x)] <- mu
    ##     (diag(nrow(A0))-t(A0))%*%nu
    ##     meanstart <- solve(diag(nrow(A0))+t(A0))%*%nu
    ##     meanstart <- meanstart[which(is.na(x$mean))]
    ##     if (length(meanstart)>0)
    ##         pp[seq(length(meanstart))] <- meanstart
    ## }
    names(pp) <- coef(x, silent=TRUE, fixed=FALSE, mean=TRUE)[seq_len(length(pp))]
    pp[!is.finite(pp) | is.nan(pp) | is.na(pp)] <- 0.01
    return(pp)
}

###}}} startvalues0

###{{{ startvalues

## McDonald & Hartmann, 1992
##' @export
startvalues <-
function(x, S, mu=NULL, debug=FALSE, silent=FALSE, tol=1e-6, delta=1e-6,...) {
  ## As proposed by McDonald & Hartmann, 1992.
  ## Implementation based on John Fox's implementation in the 'sem' R-package
  S <- reorderdata.lvm(x,S)
  if (nrow(S)!=length(manifest(x))) stop("Number of observed variables in data and models does not agree.")
  J <- index(x)$J ## Manifest selection
  P0 <- index(x)$P0 ## covariance 'adjacency'
  A <- t(index(x)$M) ## Adjacency matrix
  n <- nrow(S) ## Number of manifest variables
  m <- nrow(A) ## Number of variables
  A0 <- t(index(x)$M0) ## Adjacency matrix (without fixed parameters)
  obs.idx <- as.vector(J%*%(seq_len(m)));  latent.idx <- setdiff(seq_len(m), obs.idx)
  s <- sqrt(diag(S))
  suppressWarnings(R <- (cov2cor(S))) ## S/outer(s,s)
  C <- P0
  Debug(list("obs.idx", obs.idx), debug)
  C[obs.idx,obs.idx] <- R
  ## Estimates of covariance between latent and manifest variables
  Debug((C), debug)
  for (i in latent.idx) {
    inRelation <- A[obs.idx,i]==1
    for (j in seq_len(length(obs.idx))) {
      Debug((j), debug)
      C[obs.idx[j],i] <- C[i,obs.idx[j]] <- if (any(inRelation)) {
        numerator <- sum(R[j, which(inRelation)])
        denominator <- sqrt(sum(R[which(inRelation), which(inRelation)]))
        numerator/denominator ## as proposed by McDonald & Hartmann
      } else {
        runif(1, .3, .5) ## No arrows => small random covariance
      }
    }
  }
  ## Estimates of covariance between latent variables
  for (i in latent.idx) {
    for (j in latent.idx) {
      C[i,j] <- C[j,i] <-
        if (i==j) {
          1
        } else {
          inRelation.i <- A[obs.idx, i]==1
          inRelation.j <- A[obs.idx, j]==1
          if ((any(inRelation.i)) | (any(inRelation.j))) {
            numerator <- sum(R[which(inRelation.i), which(inRelation.j)])
            denominator <- sqrt( sum(R[which(inRelation.i), which(inRelation.i)])
                                * sum(R[which(inRelation.j), which(inRelation.j)]))
            numerator/(denominator+0.01) ## Avoid division by zero
          } else {
            runif(1, .3, .5)
          }
        }
    }
  }
  if (debug) {
    print("C="); print(C);
  }
  Ahat <- matrix(0,m,m)
  C[is.nan(C)] <- 0
  for (j in seq_len(m)) { ## OLS-estimates
    relation <- A[j,]==1
    if (!any(relation)) next
    Ahat[j, relation] <- tryCatch(Inverse(C[relation,relation] + diag(nrow=sum(relation))*delta,tol=1e-3) %*% C[relation,j], error=function(...) 0)
  }
  Ahat[obs.idx,] <- Ahat[obs.idx,]*matrix(s, n, m)
  Ahat[,obs.idx] <- Ahat[,obs.idx]/matrix(s, m, n, byrow=TRUE)
  Chat <- C
  Chat[obs.idx,] <- Chat[obs.idx,]*matrix(s,n,m)  ##
  Chat[,obs.idx] <- Chat[,obs.idx]*matrix(s,m,n,byrow=TRUE)  ##
  Phat <- (diag(m)-Ahat)%*%Chat%*%t(diag(m)-Ahat)
  ##diag(Phat) <- abs(diag(Phat))
  ## Guarantee PD-matrix:
  Phat[is.nan(Phat) | is.na(Phat)] <- 0
  diag(Phat)[diag(Phat)==0] <- 1
  eig <- eigen(Phat)
  L <- abs(eig$values); L[L<1e-3] <- 1e-3
  Phat <- eig$vectors%*%diag(L,ncol=ncol(eig$vectors))%*%t(eig$vectors)

  Debug(list("start=",start), debug)
  start <- pars(x, A=t(Ahat*A0), P=(Phat*P0))
  names(start) <- coef(x, silent=TRUE, fixed=FALSE, mean=FALSE)[seq_len(length(start))]
  res <- startmean(x,start,mu)
  res[!is.finite(res) | is.nan(res) | is.na(res)] <- 1
  res
}

###}}} startvalues
