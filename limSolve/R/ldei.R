
##==============================================================================
## ldei Solves underdetermined problem,
## Least Distance with equalities and inequalities
##==============================================================================

ldei <- function(E, F, G=NULL, H=NULL,
  tol=sqrt(.Machine$double.eps), verbose=TRUE) {

## input consistency
  if (! is.matrix(E) & ! is.null(E))
    E <- t(as.matrix(E))
  if (! is.matrix(G) & ! is.null(G))
    G <- t(as.matrix(G))

##------------------------------------------------------------------------
## 0. Setup problem
##------------------------------------------------------------------------
  if (is.null(tol))
    tol <- sqrt(.Machine$double.eps)

  ## Problem dimension
  Neq    <- nrow(E)   # number of equations
  Nx     <- ncol(E)   # number of unknowns
  Nin    <- nrow(G)   # number of inequalities

  ## copies
  EE <- E
  GG <- G

  ## Consistency of input
  if (!is.null(G)) {
    if (ncol(G)   != Nx)
      stop("cannot solve least distance problem - E and G not compatible")
    if (length(H) != Nin)
      stop("cannot solve least distance problem - G and H not compatible")
  }

  if (length(F) != Neq)
    stop("cannot solve least distance problem - E and F not compatible")

##------------------------------------------------------------------------
## 1. Decompose matrix by singular value decomposition
##------------------------------------------------------------------------
  S          <- svd(E,nrow(E),ncol(E))

  ## number of 'solvable unknowns' - the rank of the problem: 
  solvable   <- sum(S$d > tol * S$d[1])
  if (solvable > Nx)
    stop("cannot solve problem by LDEI - overdetermined system")

  ## number of 'unsolvable' unknowns
  unsolvable <- Nx-solvable

##------------------------------------------------------------------------
## 2. Backsubstitution
##------------------------------------------------------------------------
  Xb         <- S$v[,1:solvable] %*% (t(S$u[,1:solvable])/S$d[1:solvable])%*%F

  ## Check if solution is correct
  CC         <- E %*%Xb - F
  if (any(abs(CC)>tol))
    stop("cannot solve problem by LDEI - equations incompatible")

  ##------------------------------------------------------------------------
  ## 3. Constrained solution
  ##------------------------------------------------------------------------
  ## Check if inequalities are already satisfied

  ifelse (!is.null(G), CC        <- G%*%Xb-H, CC<-0)

  if (all(CC > -tol)) {
     X    <- Xb
     IsError <-FALSE
  } else   {

   ## Number of unknowns to be solved by the least distance programming
    LDPNx   <-unsolvable

   ## Construct an orthogonal basis of V2 = (I-V*VT)*Random
    rnd     <- matrix(data=runif(Nx*unsolvable),nrow=Nx,ncol=unsolvable)
    V2      <- diag(1,Nx) - S$v[,1:solvable]%*%t(S$v[,1:solvable])
    V2      <- V2 %*% rnd

   ## Orthogonal basis, constructed by Singular Value Decomposition on V2.
    Vort    <- svd(V2,nu=0,nv=unsolvable)  ## only right-hand side needed
    ortho   <- V2 %*% Vort$v[,1:unsolvable]
    for (i in 1:Nx)
      ortho[i,]<- ortho[i,]/Vort$d[1:unsolvable]

   ## Least distance programming matrices LDPG = G*ortho, LDPH = H - G*Xb
    LDPG    <- G %*% ortho
    LDPH    <- H - G %*% Xb

   ## call the DLL with the least distance programming routine ldp
    IsError <- FALSE
    NW      <- (LDPNx+1)*(Nin+2) +2*Nin

    storage.mode(LDPG) <- storage.mode(LDPH) <- "double"
    
    sol  <-.Fortran("ldp",G=LDPG,H=LDPH,
             NUnknowns=as.integer(LDPNx),NConstraints=as.integer(Nin),
             NW=as.integer(NW),X=as.vector(rep(0,LDPNx)),XNorm=0.,
             W=as.double(rep(0.,NW)),xIndex=as.integer(rep(0,Nin)),
             Mode=as.integer(0),
             verbose=as.logical(verbose),IsError=as.logical(IsError))

    IsError<-sol$IsError

   ## The solution, corrected for weights
    X    <- ortho[,1:unsolvable] %*% as.matrix(sol$X) + Xb
    X[which(abs(X)<tol)] <- 0         ## zero very tiny values
      
  } ## end CC > - tol

  ## Residual of the inequalities
  residin  <- 0
  if (!is.null(GG)) {
    ineq     <- GG %*% X - H
    residin  <- -sum(ineq[ineq<0])
  }
  ## Total residual norm
  residual <- sum(abs(EE %*% X - F))+residin

  ## The solution norm
  solution <- sum ((X)^2)
  X    <- as.vector(X)
  xnames <- colnames(E)
  if (is.null(xnames)) xnames <- colnames(G)
  names (X) <- xnames

  return(list(X=X,
              unconstrained.solution=as.vector(Xb),
              residualNorm=residual,
              solutionNorm=solution,
              IsError=IsError,
              type="ldei"))

}
