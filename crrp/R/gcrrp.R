
gcrrp <- function(time, fstatus, X, failcode=1, cencode=0, group=1:ncol(X), penalty=c("gLASSO", "gMCP", "gSCAD"), 
                  gamma=switch(penalty, SCAD=3.7, 2.7),alpha=1, lambda.min=0.001, nlambda=50, lambda, 
                  eps=.001, max.iter=1000, weighted=FALSE) {
  
  norm <- function(beta){
    temp <- sqrt(sum(beta^2))
    return(temp)
  }

  ginv = function(X, tol = sqrt(.Machine$double.eps)){
    s = svd(X)
    nz = s$d > tol * s$d[1]
    if(any(nz)) s$v[,nz] %*% (t(s$u[,nz])/s$d[nz])
    else X*0
  }
  
  dgpenalty <- function(beta, K, gamma=switch(penalty, SCAD=3.7, 2.7), lambda, penalty=c("gLASSO", "gSCAD", "gMCP"),
                       penalty.factor){
    penalty <- match.arg(penalty)
    m <- penalty.factor
    ans <- rep(NA, length(beta))
    lambda <- lambda*sqrt(K)*m
    if (penalty=="gLASSO"){
      ans[norm(beta) > 0] = (norm(beta[norm(beta) > 0]))^{-1}
      ans[norm(beta)==0]=0
      return(ans*lambda)
      
    }
    else if (penalty=="gSCAD"){
      ans[norm(beta) > 0]=(1*(norm(beta[abs(beta) > 0])<=lambda)*lambda+
                             pmax(gamma*lambda-norm(beta[abs(beta) > 0]), 0)/(gamma-1)*(norm(beta[abs(beta) > 0])>lambda))/norm(beta[abs(beta) > 0])
      ans[beta==0]=0
      return(ans)
    }
    else if (penalty=="gMCP"){
      ans[norm(beta)>0]=1*(norm(beta[norm(beta) > 0])<=gamma*lambda)*(lambda-norm(beta[norm(beta) > 0])/gamma)/norm(beta[norm(beta) > 0])
      ans[norm(beta)==0]=0
      return(ans)
      
    }
  }
  
  
  orthogonalize <- function(X, group) {
    n <- nrow(X)
    J <- max(group)
    T <- vector("list", J)
    XX <- matrix(0, nrow=nrow(X), ncol=ncol(X))
    for (j in seq_along(numeric(J))) {
      ind <- which(group==j)
      if (length(ind)==0) next
      SVD <- svd(X[, ind, drop=FALSE], nu=0)
      r <- which(SVD$d > 1e-10)
      T[[j]] <- sweep(SVD$v[,r,drop=FALSE], 2, sqrt(n)/SVD$d[r], "*")
      XX[,ind[r]] <- X[,ind]%*%T[[j]]
    }
    nz <- !apply(XX==0,2,all)
    XX <- XX[, nz, drop=FALSE]
    attr(XX, "T") <- T
    attr(XX, "group") <- group[nz]
    XX
  }
  unorthogonalize <- function(b, XX, group) {
    ind <- !sapply(attr(XX, "T"), is.null)
    T <- bdiag(attr(XX, "T")[ind])
    as.matrix(T %*% b)
  }
  
  ## Error checking
  if (gamma <= 1 & penalty=="gMCP") stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty=="gSCAD") stop("gamma must be greater than 2 for the SCAD penalty")
  if (nlambda < 2) stop("nlambda must be at least 2")
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead")
  
  
  ##sort time
  n <- length(time)
  p <- ncol(X)
  d <- data.frame(time=time, fstatus=fstatus)
  if (!missing(X)) d$X <- X
  d <- d[order(d$time),]
  ftime <- d$time
  cenind <- ifelse(d$fstatus==cencode,1,0)
  fstatus <- ifelse(d$fstatus==failcode, 1,2*(1-cenind))
  X <- d$X
  u <- do.call('survfit',list(formula=Surv(ftime,cenind)~1,data=data.frame(ftime,cenind)))  
  #### uuu is weight function
  u <- approx(c(0,u$time,max(u$time)*(1+10*.Machine$double.eps)),c(1,u$surv,0),
              xout=ftime*(1-100*.Machine$double.eps),method='constant',f=0,rule=2)
  uuu <- u$y
  #end of data preparation##################################################
 
  penalty.factor=rep(1,max(group))
  ## Reorder groups, if necessary
  xnames <- if (is.null(colnames(X))) paste("V",1:ncol(X),sep="") else colnames(X)
  colnames(X) <- xnames
  if (any(order(group) != 1:length(group)) | !is.numeric(group)) {
    reorder.groups <- TRUE
    g <- as.numeric(group)
    g.ord <- order(g)
    g.ord.inv <- match(1:length(g), g.ord)
    g <- g[g.ord]
    X <- X[,g.ord]
    penalty.factor <- penalty.factor(g.ord)
  } else {
    reorder.groups <- FALSE
    g <- group
    J <- max(g)
  }
  
  ## Set up XX, yy, lambda
  std <- .Call("standardize", X, PACKAGE="crrp")
  XX <- std[[1]]
  center <- std[[2]]
  scale <- std[[3]]
  nz <- which(scale > 1e-6)
  zg <- setdiff(unique(g), unique(g[nz]))
  XX <- XX[ ,nz, drop=FALSE]
  g <- g[nz]
  XX <- orthogonalize(XX, g)
  g <- attr(XX, "group")
  K <- as.numeric(table(g))  
  
  if (weighted) {
    init <- crr(time, fstatus, XX)$coef
    size <- c(0, cumsum(K))
    weight <- NULL
    for (jj in 1: (length(size)-1)){
      weight <- c(weight, sqrt(size[jj+1]-size[jj])/norm(init[(size[jj]+1):size[jj+1]]))
    }   
    penalty.factor <- weight
  }
  
  if (length(zg)) {
    J  <- J - length(zg)
    penalty.factor <- penalty.factor[-zg]
  }
  
  if (missing(lambda)) {
    eta0 <- rep(0,n)  
    sw <- .C("scorehessian", as.double(ftime),as.integer(fstatus), as.double(XX),
             as.integer(p), as.integer(n),  as.double(uuu), as.double(eta0), 
             double(n), double(n), double(1), PACKAGE="crrp")
    score0 <- sw[[8]]
    w0 <- sw[[9]]
    r0 <- ifelse(w0==0, 0, score0/w0)
    z <- eta0+r0
    l.max <- max(t(w0*z)%*%XX)/n
    l.min <- lambda.min
    lambda=exp(seq(log(l.max),log(l.min*l.max),len=nlambda)) 
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }
  
  ## Fit

  K1 <- as.integer(c(0, cumsum(K)))
  
  res <- .Call("gcdfit_psh", XX, as.numeric(ftime), as.integer(fstatus), uuu, K1, penalty,
               lambda, eps, as.integer(max.iter), gamma, penalty.factor, alpha, PACKAGE="crrp")
  b <- matrix(res[[1]], p, (nlambda))
  dev <- res[[2]]
  iter <- res[[3]]  
  score <- matrix(res[[5]], n, nlambda)
  hessian <- matrix(res[[6]], n, nlambda)
  
  #Selection criterion for tuning parameter

  bic <- gcv <- ps <- rep(0, nlambda)
  ps_sub <- rep(0, length(K1)-1)

  for ( l in 1:nlambda){
    zwz <- t(XX)%*%diag(hessian[,l])%*%XX
    pp <- NULL
    for ( k in 1 : (length(K1)-1)){
    pp <- c(pp,dgpenalty(beta=b[(K1[k]+1):K1[k+1],l ],K[k], lambda=lambda[l], penalty=penalty, penalty.factor=penalty.factor[k]))
    }   
    
    ps[l] <- sum(diag(XX%*%ginv(zwz+n*diag(pp))%*%t(XX)%*%diag(hessian[,l])))-sum(b[,l] == 0)      
    ll <-dev[l]/-2
    bic[l] <- -2*ll+sum(1-(b[,l ]==0))*log(n)
    gcv[l] = -ll/(n*(1-ps[l]/n)^2)
  }

 ## Unstandardize
 b <- unorthogonalize(b, XX, g)
 b <- b/scale[nz]
 beta <- matrix(0, nrow=ncol(X), ncol=length(lambda))
 if (reorder.groups) {
   beta <- b[g.ord.inv,]
 } else {
   beta<- b
 }
 
 ## Names
 
 dimnames(beta) <- list(xnames, round(lambda,digits=4))
 
  val <- structure(list(beta=beta,
                        iter=iter,
                        group=group,
                        lambda=lambda,
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        loglik = ll,
                        GCV = gcv,
                        BIC = bic),
                   class = "gcrrp")
  val
}
