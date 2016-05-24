testConstraints <- function(model, qhat, uhat, constraints, method=c("D1","D2"), df.com=NULL){
# test constraints with multiply imputed data
  
  if(missing(model)==(missing(qhat)|missing(uhat))) stop("Either 'model' or both 'qhat' and 'uhat' must be supplied.")

  # match methods
  if(!missing(method) & length(method)>1) stop("Only one of 'D1' or 'D2' may be supplied as 'method'.")
  method <- match.arg(method)

  cons <- gsub("\\(Intercept\\)","Intercept",constraints)

  # warnings for ignored arguments
  if(!is.null(df.com) & method=="D2") warning("Setting complete-data degrees of freedom is only available for 'D1' and will be ignored with 'D2'.")

  if(!missing(qhat)){

    coef.method <- "default"
    if(missing(uhat)) stop("Either 'model' or both 'qhat' and 'uhat' must be supplied.")

    if(length(dim(qhat))==3) qhat <- apply(qhat,3,identity)

    if(is.list(qhat)){
      Qhat <- sapply(qhat, identity)
      Uhat <- sapply(uhat, identity, simplify="array")
    }else{
      Qhat <- qhat
      Uhat <- uhat
    }
    if(is.null(dim(Qhat))){ 
      dim(Qhat) <- c(1,length(qhat))
      nms <- if(is.list(qhat)) names(qhat[[1]]) else if(is.matrix(qhat)) rownames(qhat) else names(qhat)[1]
      dimnames(Qhat) <- list(nms, NULL)
    }
    if(is.null(dim(Uhat))) dim(Uhat) <- c(1,1,length(qhat))
    if(is.null(rownames(Qhat))){
      nms <- !sapply(dimnames(Uhat),is.null)
      nms <- dimnames(Uhat)[[min(which(nms))]]
      rownames(Qhat) <- nms
    }
    m <- dim(Qhat)[2]
    p <- dim(Qhat)[1]
    k <- q <- length(cons)

  }

  if(!missing(model)){

    if(!"list"%in%class(model)) stop("The 'model' argument must be a list of fitted statistical models.")

    # ***
    # select extraction methods
    cls <- class(model[[1]])

    # default method
    coef.method <- "default"

    # merMod (lme4)
    if(any(grepl("merMod",cls)) & coef.method=="default"){
      if(!requireNamespace("lme4", quietly=TRUE)) stop("The 'lme4' package must be installed to handle 'merMod' class objects.")
      coef.method <- "lmer"
    }
  
    # lme (nlme)
    if(any(grepl("^.?lme$",cls)) & coef.method=="default"){
      if(!requireNamespace("nlme", quietly=TRUE)) stop("The 'nlme' package must be installed to handle '(n)lme' class objects.")
      coef.method <- "nlme"
    }

    fe <- switch(coef.method,
      lmer=.getCOEF.lmer(model),
      nlme=.getCOEF.nlme(model),
      default=.getCOEF.default(model)
    )

    m <- length(model)
    Qhat <- fe$Qhat
    Uhat <- fe$Uhat
    if(is.null(dim(Qhat))) dim(Qhat) <- c(1,m)
    if(is.null(dim(Uhat))) dim(Uhat) <- c(1,1,m)
    p <- dim(Qhat)[1]
    k <- q <- length(cons)

  }

  newQhat <- array(NA, dim=c(q,m))
  newUhat <- array(NA, dim=c(q,q,m))

  # *** delta method
  for(ii in 1:m){

    theta <- Qhat[,ii]
    Sigma <- Uhat[,,ii]
    names(theta) <- gsub("\\(Intercept\\)","Intercept",names(theta))

    g <- parse(text=cons)
    env.g <- new.env()
    for(pp in 1:p) assign(names(theta)[pp],theta[pp],pos=env.g)

    # new parameter estimates
    newtheta <- numeric(length(g))
    for(qq in 1:q) newtheta[qq] <- eval(g[qq],envir=env.g)

    # derivative, new covariance matrix
    gdash.theta <- matrix(NA,q,p)
    #for(qq in 1:q) for(pp in 1:p) gdash.theta[qq,pp] <- eval( D(g[qq],names(theta)[pp]), envir=env.g )
    for(qq in 1:q){
      tmp <- numericDeriv(g[[qq]],names(theta),env.g)
      gdash.theta[qq,] <- attr(tmp,"gradient")
    }
    newSigma <- gdash.theta %*% Sigma %*% t(gdash.theta)

    newQhat[,ii] <- newtheta
    newUhat[,,ii] <- newSigma

  }
 
  # *** aggregation
  if(method=="D1"){

    Qbar <- apply(newQhat,1,mean)
    Ubar <- apply(newUhat,c(1,2),mean)
    B <- cov(t(newQhat))
    r <- (1+m^(-1))*sum(diag(B%*%solve(Ubar)))/k
    Ttilde <- (1 + r)*Ubar
    
    # D1 (Li, Raghunathan and Rubin, 1991)
    val <- t(Qbar) %*% solve(Ttilde) %*% Qbar / k
    t <- k*(m-1)
  
    if(!is.null(df.com)){
      a <- r*t/(t-2)
      vstar <- ( (df.com+1) / (df.com+3) ) * df.com
      v <- 4 + ( (vstar-4*(1+a))^(-1) + (t-4)^(-1) * ((a^2*(vstar-2*(1+a))) / 
           ((1+a)^2*(vstar-4*(1+a)))) )^(-1)
    } else {
      if (t>4){ 
        v <- 4 + (t-4) * (1 + (1 - 2*t^(-1)) * (r^(-1)))^2
      }else{
        v <- t * (1 + k^(-1)) * ((1 + r^(-1))^2) / 2
      }
    }
    p <- 1-pf(val, k, v)

  }

  if(method=="D2"){

    Qbar <- Ttilde <- NULL
    dW <- sapply(1:m, function(z) t(newQhat[,z]) %*% solve(newUhat[,,z]) %*% newQhat[,z])

    # D2 (Li, Meng et al., 1991)
    dWbar <- mean(dW)
    r <- (1+m^(-1)) * var(sqrt(dW))
    val <- (dWbar/k - (m+1)/(m-1) * r) / (1+r)

    v <- k^(-3/m) * (m-1) * (1+r^(-1))^2
    p <- 1-pf(val, k, v)

  }

  out <- matrix(c(val,k,v,p,r),ncol=5)
  colnames(out) <- c("F.value","df1","df2","p.value","RIV")
    
  out <- list(
    call=match.call(),
    constraints=cons,
    test=out,
    Qbar=Qbar,
    T=Ttilde,
    m=m,
    adj.df=!is.null(df.com),
    df.com=df.com,
    method=method
  )

  class(out) <- "mitml.testConstraints"
  out

}

