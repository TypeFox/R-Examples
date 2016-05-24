icaimax <-
  function(X,nc,center=TRUE,maxit=100,tol=1e-6,Rmat=diag(nc),
           alg=c("newton","gradient"),fun=c("tanh","log","ext"),
           signs=rep(1,nc),signswitch=TRUE,rate=1,rateanneal=NULL){
    ###### ICA via (Fast and Robust) Information-Maximization
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: August 23, 2015
    
    ### initial checks
    X <- as.matrix(X)
    nobs <- nrow(X)
    nvar <- ncol(X)
    nc <- as.integer(nc[1])
    if(nc<1){ stop("Must set nc>=1 component.") }
    maxit <- as.integer(maxit[1])
    if(maxit<1){ stop("Must set maxit>=1 iteration.") }
    tol <- tol[1]
    if(tol<=0){ stop("Must set ctol>0.") }
    if(nc>min(nobs,nvar)){ stop("Too many components. Set nc<=min(dim(X)).") }
    if(nrow(Rmat)!=nc | ncol(Rmat)!=nc){ stop("Input 'Rmat' must be nc-by-nc rotation matrix.") }
    fun <- fun[1]
    if(fun=="ext"){
      signs <- sign(signs)
      if(length(signs)!=nc){ stop("Input 'signs' must be have length equal to 'nc' input.") }
    } else {
      signs <- NA
      signswitch <- FALSE
    }
    alg <- alg[1]
    if(alg=="gradient"){
      rate <- rate[1]
      if(rate<=0){ stop("Must set 'rate' greater than 0.") }
      if(!is.null(rateanneal[1])){
        if(length(rateanneal)!=2L){ stop("Input 'rateanneal' should be two-element vector.") }
        if(rateanneal[1]<=0 || rateanneal[1]>=90){ stop("Input 'rateanneal[1]' should be in range (0,90).") }
        if(rateanneal[2]<=0 || rateanneal[2]>1){ stop("Input 'rateanneal[2]' should be in range (0,1].") }
        ralog <- TRUE
      } else {
        ralog <- FALSE
      }
    }
    
    ### center and whiten
    if(center) X <- scale(X,scale=FALSE)
    xeig <- eigen(crossprod(X)/nobs,symmetric=TRUE)
    nze <- sum(xeig$val>xeig$val[1]*.Machine$double.eps)
    if(nze<nc){
      warning("Numerical rank of X is less than requested number of components (nc).\n  Number of components has been redefined as the numerical rank of X.")
      nc <- nze
      Rmat <- diag(nc)
    }
    Dmat <- sdiag(sqrt(xeig$val[1:nc]))
    Mprt <- tcrossprod(Dmat,xeig$vec[,1:nc])
    diag(Dmat) <- 1/diag(Dmat)
    Pmat <- xeig$vec[,1:nc]%*%Dmat
    Xw <- X%*%Pmat   # whitened data
    
    ### check if nc=1
    if(nc==1L){
      return(list(S=Xw,M=Mprt,W=t(Pmat),Y=Xw,Q=t(Pmat),R=matrix(1),
                  vafs=(sum(Mprt^2)*nobs)/sum(X^2),iter=NA,
                  alg=alg,fun=fun,signs=signs,rate=rate))
    }
    
    ### which nonlinearity
    if(fun=="log"){
      fun1d <- function(x,sgn=1){2/(1+exp(-x))-1}
      fun2d <- function(x,sgn=1){1/(cosh(x)+1)}
    } else if(fun=="ext"){
      fun1d <- function(x,sgn=1){x+tanh(x)%*%sdiag(sgn)}
      fun2d <- function(x,sgn=1){1+(1-tanh(x)^2)%*%sdiag(sgn)}
    } else {
      fun1d <- function(x,sgn=1){tanh(x)}
      fun2d <- function(x,sgn=1){1-tanh(x)^2}
    }
    
    ### which algorithm
    if(alg=="gradient"){
      
      # gradient descent
      iter <- 0
      vtol <- 1
      while(vtol>tol && iter<maxit){
        # update all components
        smat <- Xw%*%Rmat
        if(signswitch){ signs <- sign(colMeans((cosh(smat)^-2)-tanh(smat)*smat)) }
        rnew <- Rmat - rate*crossprod(Xw/nobs,fun1d(smat,signs))
        # orthgonalize
        rsvd <- svd(rnew)
        rnew <- tcrossprod(rsvd$u,rsvd$v)
        # check for convergence
        vtol <- 1 - min(abs(colSums(Rmat*rnew)))
        iter <- iter + 1
        Rmat <- rnew
        if(ralog && ((acos(1-vtol)*180/pi)<rateanneal[1])){ rate <- rate*rateanneal[2] }
      } # end while(vtol>tol && iter<maxit)
      
    } else{
      
      # Newton iteration
      iter <- 0
      vtol <- 1
      while(vtol>tol && iter<maxit){
        # update all components
        smat <- Xw%*%Rmat
        if(signswitch){ signs <- sign(colMeans((cosh(smat)^-2)-tanh(smat)*smat)) }
        Hmat <- matrix(colMeans(fun2d(smat,signs)),nc,nc,byrow=TRUE)
        rnew <- Rmat - crossprod(Xw/nobs,fun1d(smat,signs))/Hmat
        # orthgonalize
        rsvd <- svd(rnew)
        rnew <- tcrossprod(rsvd$u,rsvd$v)
        # check for convergence
        vtol <- 1 - min(abs(colSums(Rmat*rnew)))
        iter <- iter+1
        Rmat <- rnew
      } # end while(vtol>tol && iter<maxit)
      
    } # end if(alg=="gradient")
    
    ### sort according to vafs
    M <- crossprod(Rmat,Mprt)
    vafs <- rowSums(M^2)
    ix <- sort(vafs,decreasing=TRUE,index.return=TRUE)$ix
    M <- M[ix,]
    Rmat <- Rmat[,ix]
    vafs <- (vafs[ix]*nobs)/sum(X^2)
    
    return(list(S=Xw%*%Rmat,M=t(M),W=t(Pmat%*%Rmat),Y=Xw,
                Q=t(Pmat),R=Rmat,vafs=vafs,iter=iter,
                alg=alg,fun=fun,signs=signs,rate=rate))
    
  }