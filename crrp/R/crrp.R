
crrp <- function(time, fstatus, X, failcode=1, cencode=0,
                   penalty=c("MCP", "SCAD", "LASSO"), gamma=switch(penalty, SCAD=3.7, 2.7), 
                   alpha=1, lambda.min=.001, nlambda=50, lambda, eps=.001, 
                   max.iter=1000, penalty.factor=rep(1, ncol(X)), weighted=FALSE){
 
  dpenalty <- function(beta, gamma=switch(penalty, SCAD=3.7, 2.7), lambda, penalty=c("LASSO", "SCAD", "MCP"),
                       penalty.factor=rep(1, length(beta))){
    penalty <- match.arg(penalty)
    m <- penalty.factor
    ans <- rep(NA, length(beta))
    if (penalty=="LASSO"){
      ans[abs(beta) > 0] = (abs(beta[abs(beta) > 0])*m[abs(beta)>0])^{-1}
      ans[beta==0]=0
      w <- diag(as.vector(ans))
      ww <- w*lambda
      return(ww)
    }
    else if (penalty=="SCAD"){
      ans[abs(beta) > 0]=(lambda*(abs(beta[abs(beta) > 0])<=lambda)+
                            pmax(gamma*lambda-abs(beta[abs(beta) > 0]), 0)/(gamma-1)*(abs(beta[abs(beta) > 0])>lambda))/abs(beta[abs(beta) > 0])
      ans[beta==0]=0
      return(diag(as.vector(ans)))
    }
    else if (penalty=="MCP"){
      ans[abs(beta)>0]=1*(abs(beta[abs(beta) > 0])<=gamma*lambda)*(lambda-abs(beta[abs(beta) > 0])/gamma)/abs(beta[abs(beta) > 0])
      ans[beta==0]=0
      return(diag(as.vector(ans)))
      
    }
  }
  
  ginv = function(X, tol = sqrt(.Machine$double.eps)){
    s = svd(X)
    nz = s$d > tol * s$d[1]
    if(any(nz)) s$v[,nz] %*% (t(s$u[,nz])/s$d[nz])
    else X*0
  }
  
  ##sort time
    n <- length(time)
    p <- ncol(X)
    d <- data.frame(time=time, fstatus=fstatus)
    if (!missing(X)) d$X <- as.matrix(X)
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
  #Output time, fstatus, uuu, X
  
  ## Error checking
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
  if (nlambda < 2) stop("nlambda must be at least 2")
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead")
  if (length(penalty.factor)!=ncol(X)) stop("penalty.factor does not match up with X")

  ## Set up XX, yy, lambda
  std <- .Call("standardize", X, PACKAGE="crrp")
  XX <- std[[1]]
  center <- std[[2]]
  scale <- std[[3]]
  
  if (weighted) penalty.factor <- penalty.factor*scale
  
  nz <- which(scale > 1e-6)
  if (length(nz) != ncol(XX)) XX <- XX[ ,nz, drop=FALSE]
  penalty.factor <- penalty.factor[nz]
  
  # initial value
  #set value of lambda
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

  res <- .Call("cdfit_psh", XX, as.numeric(ftime), as.integer(fstatus), uuu, penalty, lambda, 
               eps, as.integer(max.iter), gamma, penalty.factor, alpha, PACKAGE="crrp")
  b <- matrix(res[[1]], p, (nlambda))
  dev <- res[[2]]
  iter <- res[[3]]  
  score <- matrix(res[[5]], n, nlambda)
  hessian <- matrix(res[[6]], n, nlambda)  
  ## Unstandardize
  beta <- b/scale  
  ## Names
  if (is.null(colnames(X))) varnames <- paste("V",1:ncol(X),sep="")
  else varnames <- colnames(X)
  colnames(beta) <-  round(lambda,digits=4)
  rownames(beta) <- varnames
  
  #calculate GCV, BIC, and standard error
  sd <- matrix(0, p, nlambda)
  bic <- gcv <- ps <- rep(0, nlambda)
  for ( l in 1:nlambda){
  zwz <- t(XX)%*%diag(hessian[,l])%*%XX
  pp <- dpenalty(beta=b[,l ],lambda=lambda[l], penalty=penalty, penalty.factor=penalty.factor)
  sd[,l ] <- sqrt(diag(ginv(zwz+n*pp)%*%zwz%*%ginv(zwz+n*pp)))/scale
  
  ll <-dev[l]/-2
  bic[l] <- -2*ll+sum(1-(b[,l ]==0))*log(n)
  
  ps[l] <- sum(diag(XX%*%ginv(zwz+n*pp)%*%t(XX)%*%diag(hessian[,l])))-sum(b[,l ] == 0)
  gcv[l] = -ll/(n*(1-ps[l]/n)^2)
  }
  
  ## Output
  val <- structure(list(beta = beta,
                        iter = iter,
                        lambda = lambda,
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        loglik = ll,
                        GCV = gcv,
                        BIC = bic,
                        SE=sd,
                        call=sys.call()),
                   class = "crrp")

  val
  
}
  