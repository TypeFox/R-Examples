#Predict.matrix <- mgcv::Predict.matrix
#formals(Predict.matrix) <- c(formals(Predict.matrix), alist(... = ))

#### Define penalized splines with truncated lines
#### m is order of the spline, m-1 is the degree 

#p-splines with truncated polynomials
`smooth.construct.tl.smooth.spec` <-
function(object,data,knots)
{ 
  m<-object$p.order # p+1
  if (is.na(m)) m<-3
  if (m==1)  {warning("for estimation the order of the spline was set to its default 3"); m <- 3}
  m <- c(m-1,m)
  nk<-object$bs.dim  # number of interior knots
  x <- get.var(object$term,data)  # find the data
  if (!is.null(knots)) k <- get.var(object$term,knots) 
  else k<-NULL
  if (is.null(k)) 
  k <- quantile(unique(x),seq(0,1,length=(nk+2))[-c(1,(nk+2))])
  names(k) <- NULL   
  z=outer(x,k,"-")
  z=z*(z>0)
  Z=z^(m[1])
  X=rep(1,length(x)) 
  for (i in 1:(m[1])) X=cbind(X,x^i)
  object$X<-cbind(X,Z) # get model matrix
  if (!object$fixed)       
  { S<-diag(c(rep(0,m[2]),rep(1,nk)))
    object$S<-list(S)  # get penalty
     }
  object$rank<-nk  # penalty rank 
  object$null.space.dim <- m[1]  # dimension of unpenalized space  
   object$knots<-k;object$m<-m      # store p-spline specific info.
   object$df<-ncol(object$X)     # maximum DoF
 
  class(object)<-"tlspline.smooth"  # Give object a class
  object
}
 
`Predict.matrix.tlspline.smooth` <-
function(object,data,...)
# prediction method function for the p.spline smooth class
{ 
  x <- get.var(object$term,data)
  m=object$m
  k=object$knots
  z=outer(x,k,"-")
  z=z*(z>0)
  Z=z^(m[1])
  X=rep(1,length(x))   
  for (i in 1:(m[1])) X=cbind(X,x^i)
  X<-cbind(X,Z) # get model matrix
  X
}


#### Define penalized splines with cubic B-splines 
#### and the penalty as the integrated squared second derivative
#### see Wand and Ormerod (2008) - O'Sullivan splines

`smooth.construct.os.smooth.spec` <- function(object,data,knots)

{  
    m <- object$p.order
    if (length(m)<2) if (is.na(m)) m <- c(3,2)
    object$p.order <- m
    if (length(m)==1){
      if (((object$p.order+1)/2)%%1 != 0) {warning("If only degree of B-spline basis is given, degree must be chosen such that q=(p+1)/2 is an integer. Set to its default p=3, q=2."); m =c(3,2)}
      else m <- c(object$p.order, (object$p.order+1)/2)
    }

    p=m[1]
    q=m[2]

    if (object$bs.dim < 0) object$bs.dim <- max(10, m[1] + 1)
    nk <- object$bs.dim
    if (nk <= 0) stop("basis dimension too small for b-spline order")

    x <- data[[object$term]]

#if (any(round(diff(knots[[object$term]]),6)!=round(diff(knots[[object$term]])[1],6))) warning("Forcing equidistant knots.")

    knots= seq(min(x),max(x),length=nk+2)[-c(1,nk+2)]
    names(knots) <- NULL

    object$X <-bs(x,knots=knots,degree=m[1],Boundary.knots=c(min(x),max(x)),intercept=T)
    object$X=t(apply(object$X,1,function(x) x-colSums(object$X)/length(data[[object$term]]))) #      C=(diag(rep(1,n))-1/n)%*%C

    d=diag(ncol(object$X))
    for(i in 1:q){
      d=diff(d)
      allKnots_p <- c(rep(min(x),p+1-i),knots,rep(max(x),p+1-i))
      weights=matrix(rep(1/(allKnots_p[-(1:(p+1-i))]-allKnots_p[-((nk+2*p-p-i+2):(nk+2*p))]),each=ncol(d)),ncol(d),nrow(d))
      d=d*(p+1-i)*t(weights)
    }

    #calculate the matrix R (integral over product of two B-spline bases)
      allKnots <- c(rep(min(x),(p-q)+1),knots,rep(max(x),(p-q)+1))
      R_int= matrix(0,nk+p+1-q,nk+p+1-q)
      for (i in 1:(nk+p-q+1)){
        for (j in i:(min(i+p-q,nk+p-q+1))){
          R <- function(x){
            Nq <- spline.des(allKnots,x,p-q+1,derivs=0*x,outer.ok=T)$design
            Nq[,i]*Nq[,j]
          }
          x1=allKnots[j]
          x2=allKnots[(i+(p-q)+1)]
          R_int[i,j] <- integrate(R,x1,x2,subdivisions=1000)$value
        }
      }
      R_int <- R_int+t(R_int)
      diag(R_int)<- diag(R_int)/2

      Dq=t(d)%*%R_int%*%d

    object$S<-list(Dq)
    object$rank <- object$bs.dim - m[2]
    object$null.space.dim <- m[2]
    object$knots <- knots
    object$m <- m

    class(object) <- "ospline.smooth"
    object
}


Predict.matrix.ospline.smooth<-function(object,data,drv=0,...){
# prediction method function for the p.spline smooth class
  X <-bs2(data[[object$term]],knots=object$knots,degree=object$m[1],Boundary.knots=c(min(data[[object$term]]),max(data[[object$term]])),intercept=T,drv=drv)
  X
}



Predict.matrix.lme <-function(object,data,drv=0,center=T,...)
# prediction method function for the p.spline smooth class
{   
  m <- object$m
  q=m[2]
  p=m[1]
  x <- data[[object$term]]
  n=length(x)
  k=object$knots
  nk <- object$bs.dim
  if (class(object)=="ospline.smooth"){
    X <- Predict.matrix.ospline.smooth(object, data,drv=drv)  # Model matrix

    d=diag(ncol(X))
    for(i in 1:q){
      d=diff(d)
      allKnots_p <- c(rep(min(x),p+1-i),k,rep(max(x),p+1-i))
      weights=matrix(rep(1/(allKnots_p[-(1:(p+1-i))]-allKnots_p[-((nk+2*p-p-i+2):(nk+2*p))]),each=ncol(d)),ncol(d),nrow(d))
      d=d*(p+1-i)*t(weights)
    }

    #calculate the matrix R (integral over product of two B-spline bases)
      allKnots <- c(rep(min(x),(p-q)+1),k,rep(max(x),(p-q)+1))
      R_int= matrix(0,nk+p+1-q,nk+p+1-q)
      for (i in 1:(nk+p-q+1)){
        for (j in i:(min(i+p-q,nk+p-q+1))){
          R <- function(x){
            Nq <- spline.des(allKnots,x,p-q+1,derivs=0*x,outer.ok=T)$design
            Nq[,i]*Nq[,j]
          }
          x1=allKnots[j]
          x2=allKnots[(i+(p-q)+1)]
          R_int[i,j] <- integrate(R,x1,x2,subdivisions=1000)$value
        }
      }
      R_int <- R_int+t(R_int)
      diag(R_int)<- diag(R_int)/2

    # calculate the Durban-decomposition for the mixed models
      Re=eigen(R_int)
      Re12=Re$vectors%*%diag(sqrt(Re$values))%*%t(Re$vectors)
      D=Re12%*%d
      DI=tcrossprod(D)  #DI=D%*%t(D)

    # centering
      if (center) X=t(apply(X,1,function(x) x-colSums(X)/n))

    Z=X%*%t(D)%*%solve(DI)
    Dq=t(d)%*%R_int%*%d
    O.e=eigen(Dq)
    null.space=(ncol(X)-q+1):(ncol(X)-1)
    U0=O.e$vectors[,null.space]
    C=X%*%U0
    dimnames(C)[[2]]=NULL

    newknots=seq(min(x),max(x),length=nk+p+1-q+2)[-c(1,nk+p+1-q+2)]
    return(list(C=C,Z=Z, knots=newknots))
  }
  else if (class(object)=="tlspline.smooth" | class(object)=='trunc.poly') {
    x<- as.vector(data[[object$term]])
    Z <- outer(x,k,"-")
    Z <- (Z*(Z>0))^m[1]
    C=rep(1,length(x))
    for (i in 1:(m[1])) C=cbind(C,x^i)
    #centering
      if(center){
      n=nrow(C)
      colSC= colSums(C)
      colSZ= colSums(Z)
      C=t(apply(C,1,function(x) x-colSC/n)) #C=(diag(rep(1,n))-1/n)%*%C
      Z=t(apply(Z,1,function(x) x-colSZ/n)) #Z=(diag(rep(1,n))-1/n)%*%Z
      }
  }
  else if (class(object)=="tps"|class(object)=='ts.smooth'){
    if (is.null(m)) m=object$p.order
    if (is.null(k)) stop("No knots given in smooth.construct.")
    x<- as.vector(data[[object$term]])
    svd.Omega = svd(abs(outer(k,k,"-"))^m[1])
    matrix.sqrt.Omega = t(svd.Omega$v %*% (t(svd.Omega$u) * sqrt(svd.Omega$d)))
    Z = t(solve(matrix.sqrt.Omega, t(abs(outer(x, k,"-")^m[1]))))
    C = cbind(rep(1,length(x)))
    for (i in 1:((m[1]-1)/2)) C=cbind(C,x^i)
    #centering
      if(center){
        n=nrow(C)
        colSC= colSums(C)
        colSZ= colSums(Z)
        C=t(apply(C,1,function(x) x-colSC/n))  #C=(diag(rep(1,n))-1/n)%*%C
        Z=t(apply(Z,1,function(x) x-colSZ/n))  #Z=(diag(rep(1,n))-1/n)%*%Z
      }
  }
  else stop ("scbM can be fitted only with os, tl or tps basis functions")
  list(C=C,Z=Z)
}



bs2=function (x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,     Boundary.knots = range(x),drv=0)
{
    nx <- names(x)
    x <- as.vector(x)
    nax <- is.na(x)
    if (nas <- any(nax))
        x <- x[!nax]
    if (!missing(Boundary.knots)) {
        Boundary.knots <- sort(Boundary.knots)
        outside <- (ol <- x < Boundary.knots[1L]) | (or <- x >
            Boundary.knots[2L])
    }
    else outside <- FALSE
    ord <- 1 + (degree <- as.integer(degree))
    if (ord <= 1)
        stop("'degree' must be integer >= 1")
    if (!missing(df) && missing(knots)) {
        nIknots <- df - ord + (1 - intercept)
        if (nIknots < 0) {
            nIknots <- 0
            warning("'df' was too small; have used  ", ord -
                (1 - intercept))
        }
        knots <- if (nIknots > 0) {
            knots <- seq.int(from = 0, to = 1, length.out = nIknots +
                2)[-c(1, nIknots + 2)]
            stats::quantile(x[!outside], knots)
        }
    }
    Aknots <- sort(c(rep(Boundary.knots, ord), knots))
    if (any(outside)) {
        warning("some 'x' values beyond boundary knots may cause ill-conditioned bases")
        derivs <- 0:degree
        scalef <- gamma(1L:ord)
        basis <- array(0, c(length(x), length(Aknots) - degree -
            1L))
        if (any(ol)) {
            k.pivot <- Boundary.knots[1L]
            xl <- cbind(1, outer(x[ol] - k.pivot, 1L:degree,
                "^"))
            tt <- spline.des(Aknots, rep(k.pivot, ord), ord,
                derivs)$design
            basis[ol, ] <- xl %*% (tt/scalef)
        }
        if (any(or)) {
            k.pivot <- Boundary.knots[2L]
            xr <- cbind(1, outer(x[or] - k.pivot, 1L:degree,
                "^"))
            tt <- spline.des(Aknots, rep(k.pivot, ord), ord,
                derivs)$design
            basis[or, ] <- xr %*% (tt/scalef)
        }
        if (any(inside <- !outside))
            basis[inside, ] <- spline.des(Aknots, x[inside],
                ord)$design
    }
    else basis <- spline.des(Aknots, x, ord,derivs=rep(drv,length(x)))$design
    if (!intercept)
        basis <- basis[, -1L, drop = FALSE]
    n.col <- ncol(basis)
    if (nas) {
        nmat <- matrix(NA, length(nax), n.col)
        nmat[!nax, ] <- basis
        basis <- nmat
    }
    dimnames(basis) <- list(nx, 1L:n.col)
    a <- list(degree = degree, knots = if (is.null(knots)) numeric(0L) else knots,
        Boundary.knots = Boundary.knots, intercept = intercept)
    attributes(basis) <- c(attributes(basis), a)
    class(basis) <- c("bs", "basis")
    basis
}


drvbasis<-function(x,degree,knots,drv=0,basis){
 if (drv>=degree) stop("WARNING: drv >= degree, derivative order must be lower than degree of the spline.")
  X <- NULL
  if (basis=="tlspline.smooth") ncol.X <- degree
  else if (basis == "tps"|basis=='ts.smooth')      ncol.X <- (degree - 1)/2
  else stop("Unsupported basis function.")

  for (pow in 1:ncol.X) {
    pow.drv <- pow - drv
    if (pow.drv >= 0) new.col <- prod(pow:(pow.drv + 1)) * (x^pow.drv)
    else new.col <- rep(0, length(x))
    X <- cbind(X, new.col)
  }
  X=cbind(0,X)

  Z <- outer(as.vector(x), knots, "-")
  if (basis=="tlspline.smooth"){
    if (degree >= drv) {
      mfac <- prod((degree - drv + 1):degree)
      Z <- mfac * (Z * (Z >    0))^(degree - drv)
    }
  }
  else if (basis=="tps"|basis=='ts.smooth'){
    if (degree >= drv) {
      mfac <- prod((degree - drv + 1):degree)
      Z <- mfac * (Z^(degree -
      drv - 1) * abs(Z))
      svd.Omega = svd(abs(outer(knots,knots,"-"))^degree)
      sqrt.Omega = t(svd.Omega$v %*% (t(svd.Omega$u) * sqrt(svd.Omega$d)))
      Z <- t(solve(sqrt.Omega, t(Z)))
    }
  }

  X <-cbind(X,Z)
  X
}


