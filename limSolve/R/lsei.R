
##==============================================================================
## lsei        : Least Squares with Equalities and Inequalities
## Solves an lsei inverse problem,
## lsei = Least Squares with Equality and Inequality constraints
## Minimise             ||A*X-B||
## subject to             E*X=F
##                        G*X>=H
## uses either the LINPACK package lsei or solve.QP from package quadprog
##==============================================================================

lsei <- function(A=NULL, B=NULL, E=NULL, F=NULL, G=NULL, H=NULL,
    Wx=NULL, Wa=NULL, type = 1, tol=sqrt(.Machine$double.eps),
    tolrank=NULL, fulloutput = FALSE, verbose=TRUE)  {

  ##------------------------------------------------------------------------
  ## Setup problem
  ##------------------------------------------------------------------------
  ## input consistency

  if (is.vector(E) & length(F)==1)
    E <- t(E)
  else if (! is.matrix(E) & ! is.null(E))
    E <- as.matrix(E)
  if (is.vector(A) & length(B)==1)
    A <- t(A)
  else if (! is.matrix(A) & ! is.null(A))
    A <- as.matrix(A)
  if (is.vector(G) & length(H)==1)
    G <- t(G)
  else if (! is.matrix(G) & ! is.null(G))
    G <- as.matrix(G)
  if (! is.matrix(F) & ! is.null(F))
    F <- as.matrix(F)
  if (! is.matrix(B) & ! is.null(B))
    B <- as.matrix(B)
  if (! is.matrix(H) & ! is.null(H))
    H <- as.matrix(H)
  if (is.null(A) && is.null (E)) {
    if(is.null(G))
      stop("cannot solve least squares problem - A, E AND G are NULL")
    A <- matrix(data=0, nrow=1,ncol=ncol(G))
    B <- 0
  }  else if (is.null(A)) {
    A <- matrix(data=E[1,], nrow=1)
    B<-F[1]
  }

  ## Problem dimension
  Neq  <- nrow(E)
  Napp <- nrow(A)
  Nx   <- ncol(A)
  Nin  <- nrow(G)
  if (is.null (Nx))
    Nx  <- ncol(E)
  if (is.null (Nx))
     Nx  <- ncol(G)
  if (is.null (Nin))
    Nin  <- 1

  ## If equalities/inequalities absent: use type=2 instead
  if (is.null (Neq)) {
    Neq <- 0
    if (verbose & type==1)
      warning("No equalities - setting type = 2")
    type = 2
    F <- NULL
  } else  {
    if (ncol(E)   != Nx)
      stop("cannot solve least squares problem - A and E not compatible")
    if (length(F) != Neq)
      stop("cannot solve least squares problem - E and F not compatible")
  }

  if (is.null(G))
    G <- matrix(data=0, nrow=1,ncol=Nx)
  if (is.null(H))
    H <- 0

  if (ncol(G)   != Nx)
    stop("cannot solve least squares problem - A and G not compatible")
  if (length(B) != Napp)
    stop("cannot solve least squares problem - A and B not compatible")
  if (length(H) != Nin)
    stop("cannot solve least squares problem - G and H not compatible")

  if (! is.null(Wa))  {
    if (length(Wa) != Napp)
      stop ("cannot solve least squares problem - Wa should have length = number of rows of A")
    A <- A*Wa
    B <- B*Wa
  }

  Tol <- tol 
  if (is.null(Tol)) Tol <- sqrt(.Machine$double.eps)

  ##------------------------------------------------------------------------
  ## Solution
  ##------------------------------------------------------------------------

  IsError <- FALSE

  if (type == 1)     {    # use LINPACKs lsei
    ineq <- Nin+Nx
    mIP  <- ineq+2*Nx+2

    ## extra options?
    lpr <- 1
    if (fulloutput)
      lpr <- lpr+3
    if (! is.null(tolrank))
      lpr <- lpr + 6
    if (! is.null(Wx))   {
      lw  <- length (Wx)
      lpr <- lpr + 2 + lw
    }
      
    ProgOpt <- rep(1.0,lpr)
    if (lpr>1)   {
      ipr <- 1
      if (fulloutput) {ProgOpt[ipr:(ipr+2)] <- c(ipr+3,1,1); ipr <- ipr+3}
      if (! is.null(tolrank))   {
        if (length(tolrank) == 1)
          tolrank <- rep(tolrank,len=2)
        ProgOpt[ipr:(ipr+5)] <- c(ipr+6,4,tolrank[1],ipr+6,5,tolrank[2])
        ipr <- ipr+6
      }
      if (! is.null(Wx))    {
        lw  <- length (Wx)
        if (lw ==1) {
          ProgOpt[ipr:(ipr+2)]<- c(ipr+3,2,1)
        } else {
          if (lw != Nx)
            stop("cannot solve least squares problem - number of weighs should be =1 or =number of unknowns")
          lw <- lw + ipr+1
          ProgOpt[ipr:lw] <- c(lw+1,3,Wx)
        }
      }
    }
    mdW <- Neq + Napp + ineq
    if (fulloutput)
      mdW <- max(mdW, Nx)
    mWS <- 2*(Neq+Nx)+max(Napp+ineq,Nx)+(ineq+2)*(Nx+7)
    storage.mode(A) <- storage.mode(B)  <- "double"
    storage.mode(E) <- storage.mode(F)  <- "double"
    storage.mode(G) <- storage.mode(H) <- "double"
    sol <-.Fortran("lsei",NUnknowns=Nx,NEquations=Neq,NConstraints=Nin,
             NApproximate=Napp,A=A,B=B,E=E,F=F,G=G,H=H,X=as.vector(rep(0,Nx)),
             mIP=as.integer(mIP),mdW=as.integer(mdW),mWS=as.integer(mWS),
             IP=as.integer(rep(0,mIP)),
             W=as.double(matrix(data=0., nrow=mdW, ncol=Nx+1)),
             WS=as.double(rep(0.,mWS)),
             lpr=as.integer(lpr),ProgOpt=as.double(ProgOpt),
             verbose=as.logical(verbose),IsError=as.logical(IsError))
    if (any(is.infinite(sol$nX)))
      sol$IsError<-TRUE
    if (fulloutput) {
      covar<-matrix(data=sol$W,nrow=mdW,ncol=Nx+1)[1:Nx,1:Nx]
      RankEq <- sol$IP[1]
      RankApp <- sol$IP[2]
    }
  } else if (type == 2)  {                 ## use R's solve.QP, package quadprog
    if (! is.null(Wx))
      stop ("cannot solve least squares problem - weights not implemented for type 2")
    if (! is.null(Wa))
      stop ("cannot solve least squares problem - weights not implemented for type 2")

    dvec  <- crossprod(A, B) # was: t(A)%*%B
    Dmat  <- crossprod(A, A) # t(A) %*% A
    diag(Dmat) <- diag(Dmat)+1e-8
    Amat  <- t(rbind(E,G))
    bvec  <- c(F,H)

    sol   <- solve.QP(Dmat ,dvec, Amat , bvec, meq=Neq)
    sol$IsError <- FALSE
    sol$X <- sol$solution

  } else
    stop("cannot solve least squares problem - type unknown")

  X <- sol$X
  X[which(abs(X)<Tol)] <- 0         # zero very tiny values

  ## Residual of the inequalities
  if (any(is.infinite(X)))  {
    residual <- Inf
    solution <- Inf
  }  else  {
    residual <- 0
    if (Nin > 0)   {
      ineq     <- G %*% X - H
      residual  <- residual -sum(ineq[ineq<0])
    }

   ## Total residual norm
   if (Neq> 0)
     residual <- residual + sum(abs(E %*% X - F))
   if (residual>Tol)
     sol$IsError<- TRUE

   ## The solution norm
     
   solution <- 0
   if (Napp > 0)
     solution <- sum ((A %*% X - B)^2)
  }
  xnames <- colnames(A)
  if (is.null(xnames))
    xnames <- colnames(E)
  if (is.null(xnames))
    xnames <- colnames(G)
  names (X) <- xnames

  res <- list(X=X,                            # vector containing the solution of the least squares problem.
              residualNorm=residual,          # scalar, the sum of residuals of equalities and violated inequalities
              solutionNorm=solution,          # scalar, the value of the minimised quadratic function at the solution
              IsError=sol$IsError,            # if an error occurred
              type="lsei")

  if (fulloutput && type == 1)  {
    res$covar<-covar
    res$RankEq <- sol$IP[1]
    res$RankApp <- sol$IP[2]
  }

  return(res)

}
