
##==============================================================================
## varranges    : Calculates ranges of inverse equations (variables)
## Given the linear constraints
##                        E*X=F
##                        G*X>=H
## and a set of variables described by the linear equations Var = EqA*X+EqB
##
## finds the minimum and maximum values of the variables
## by successively minimising and maximising each linear combination,
## using linear programming
## uses lpSolve - may fail (if frequently repeated)
##==============================================================================

varranges <- function(E=NULL, F=NULL, G=NULL, H=NULL,
    EqA, EqB=NULL, ispos=FALSE, tol=1e-8)  {

  ## input consistency
  if (! is.matrix(E) & ! is.null(E))
    E <- t(as.matrix(E))
  if (! is.matrix(G) & ! is.null(G))
    G <- t(as.matrix(G))
  if (! is.matrix(EqA) & ! is.null(EqA))
    EqA <- t(as.matrix(EqA))

  ## Dimensions of the problem
  Neq    <- nrow(E)    # number of equations
  Nx     <- ncol(E)    # number of unknowns
  Nineq  <- nrow(G)    # number of inequalities

  if (is.null(Nineq))
    Nineq <- 0
  if (is.null(Neq))
    Neq <- 0

  NVar   <- nrow(EqA)  # number of equations to minimise/maximise
  ## con: constraints ; rhs: right hand side
  ## First the equalities 

  con   <- E
  rhs   <- F
  dir   <- rep("==",Neq)
  if (Nineq > 0) {
    con   <- rbind(con,G)
    rhs   <- c(rhs,H)
    dir   <- c(dir,rep(">=",Nineq))
  }
  Range <- matrix(ncol=2,nrow=NVar,NA)

  if (ispos) {

    obj   <- vector(length = Nx)

    for (i in 1:NVar)  {
      obj        <- EqA[i,]
      lmin       <- lp("min",obj,con,dir,rhs)
      if (lmin$status == 0)
        Range[i,1] <- lmin$objval else
      if (lmin$status == 3)
        Range[i,1] <- -1e30       else
      Range[i,1] <- NA
      lmax       <- lp("max",obj,con,dir,rhs)
      if (lmax$status == 0)
        Range[i,2] <- lmax$objval else
      if (lmax$status == 3)
        Range[i,2] <- 1e30        else
      Range[i,2]  <- NA
    }
  } else {
    ## First test if problem is solvable...
    Sol <- lsei(E=E,F=F,G=G,H=H)
    if (Sol$residualNorm > tol)  {
      warning (paste("cannot proceed: problem not solvable at requested tolerance",tol))
      return(Range)
    }
    ## double the number of unknowns: x -> x1 -x2, x1>0 and x2>0
    con <- cbind(con,-1*con)
    EqA <- cbind(EqA,-1*EqA)

    for (i in 1:NVar) {
      obj <- EqA[i,]
      lmin <- lp("min", obj, con, dir, rhs)
      if(lmin$status == 0)
        Range[i, 1] <- lmin$objval else
      if(lmin$status == 3)
        Range[i, 1] <- -1e30 else
      Range[i, 1] <- NA
      lmax <- lp("max", obj, con, dir, rhs)
      if(lmax$status == 0)
        Range[i, 2] <- lmax$objval else
      if(lmax$status == 3)
        Range[i, 2] <- 1e30 else
      Range[i, 2] <- NA
    }
  }

  if (!is.null(EqB)) {
    Range[,1]<-Range[,1]-EqB
    Range[,2]<-Range[,2]-EqB
  }
  colnames(Range) <- c("min","max")
  rownames(Range) <- rownames(EqA)
  return(Range)    # a 2-column matrix with the minimum and maximum value of each equation (variable)

}
