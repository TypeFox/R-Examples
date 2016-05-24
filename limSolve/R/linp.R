
##==============================================================================
## linp        :
## Solves a linear programming problem,
## Minimise               Cost
## subject to             E*x=f
##                        G*x>=h
##
## Note: uses lp from package lpSolve
## This R-code sometimes fails and terminates R
##      for very small problems that are repeated frequently...
##==============================================================================

linp <- function(E=NULL, F=NULL, G=NULL, H=NULL, Cost,
                 ispos=TRUE, int.vec=NULL, verbose=TRUE, ...)  {

  ## input consistency
  if (! is.matrix(E) & ! is.null(E))
    E <- t(as.matrix(E))
  if (! is.matrix(G) & ! is.null(G))
    G <- t(as.matrix(G))


  ## problem dimension
  Neq  <- nrow(E)   # Number equalities
  Nx   <- ncol(E)   # Number unknowns
  Nin  <- nrow(G)   # Number inequalities

  if (is.null(Nx )) Nx  <- ncol(G)
  if (is.null(Nin)) Nin <- 0
  if (is.null(Neq)) Neq <- 0

  NumEq <- Nin+Neq  # total number of equations

  ## consistency of input
  if (!is.null(G)) {
    if (ncol(G)   != Nx)
      stop("cannot solve linear programming problem - E and G not compatible")
    if (length(H) != Nin)
      stop("cannot solve linear programming problem - G and H not compatible")
  }

  if (!is.null(E)) {
    if (length(F) != Neq)
      stop("cannot solve linear programming problem - E and F not compatible")
  }

  if (length(Cost)!= Nx)
    stop("cannot solve linear programming problem - Cost not compatible")

  IsError <- FALSE

  ## con: constraints ; rhs: right hand side

  ## the equalities:
  con   <- E
  rhs   <- F
  dir   <- rep("==",Neq)

  ## inequalities:
  if (Nin > 0)  {
    con   <- rbind(con,G)
    rhs   <- c(rhs,H)
    dir   <- c(dir,rep(">=",Nin))
  }

  if (!ispos)  {
    con  <- cbind(con, -con)
    Cost <- c(Cost, -Cost)
    if (! is.null(int.vec))
      int.vec<-c(int.vec,int.vec+Nx)
  }

  ## the solution
  sol    <- lp("min",Cost,con,dir,rhs,int.vec=int.vec,...)
  mode   <- sol$status
  ## from the original code lp_lib.h
  if (mode == -5) {
    print("unknown error")
    IsError<-TRUE
  }
  else if (mode == -4) {
    print("data ignored")
    IsError<-TRUE
  }
  else if (mode == -3) {
    print("no bfp")
    IsError<-TRUE
  }
  else if (mode == -2) {
    print("no memory")
    IsError<-TRUE
  }
  else if (mode == -1) {
    print("problem not run")
    IsError<-TRUE
  }
  else if (mode == 2) {
    print("problem infeasible")
    IsError<-TRUE
  }
  else if (mode == 3) {
    print("problem unbounded")
    IsError<-TRUE
  }
  else if (mode == 4) {
    print("problem degenerate")
    IsError<-TRUE
  }
  else if (mode == 5) {
    print("problem failed")
    IsError<-TRUE
  }
  else if (mode == 6) {
    print("user aborted")
    IsError<-TRUE
  }
  else if (mode == 7) {
    print("timed out")
    IsError<-TRUE
  }
  else if (mode == 8) {
    print("running")
    IsError<-TRUE
  }
  else if (mode == 9) {
    print("presolved")
    IsError<-TRUE
  }

  X           <- sol$solution
  if (!ispos)
    X <- X[1:Nx]-X[(Nx+1):(2*Nx)]
  solutionNorm<- sol$objval

  ## Total residual norm
  residual <- 0
  if (!is.null(E))
    residual <- sum(abs(E %*% X - F))+residual
  if (!is.null(G))   {
    ineq     <- G %*% X - H
    residual <- residual -sum(ineq[ineq<0])
  }
  xnames <- colnames(E)
  if (is.null(xnames))
    xnames <- colnames(G)
  if (is.null(xnames))
    xnames <- names(Cost[1:Nx])
  names (X) <- xnames

  return(list(X=X,                        # vector containing the solution of the linear programming problem.
              residualNorm=residual,      # scalar, the sum of residuals of equalities and violated inequalities
              solutionNorm=solutionNorm,  # scalar, the value of the minimised cost function
              IsError=IsError,            # if an error occurred
              type="linp"))

}
