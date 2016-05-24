
##==============================================================================
## xranges: estimates ranges of lsei problem
## Given the linear constraints
##                        E*X=F
##                        G*X>=H
##
## finds the minimum and maximum values of all elements of vector X
## by successively minimising and maximising each x, using linear programming
## uses lpSolve - may fail (if frequently repeated)
## unknowns can possibly be negative unless ispos=TRUE
## if all are positive, then it is solved much faster.
##==============================================================================

xranges  <-  function (E = NULL, F = NULL, G = NULL, H = NULL,
    ispos=FALSE, tol = 1e-8, central = FALSE, full=FALSE)  {

  ## input consistency
  if (!is.matrix(E) & !is.null(E))
    E <- t(as.matrix(E))
  if (!is.matrix(G) & !is.null(G))
    G <- t(as.matrix(G))
  ## Dimensions of the problem
  Neq <- nrow(E)
  Nx  <- ncol(E)
  if (is.null(Nx)) Nx <- ncol(G)
  Nineq <- nrow(G)
  if (is.null(Nineq))
    Nineq <- 0
  if (is.null(Neq))
    Neq <- 0
  Range <- matrix(ncol = 2, nrow = Nx, NA)

  ## con: constraints ; rhs: right hand side
  ## First the equalities 
  con <- E
  rhs <- F
  dir <- rep("==", Neq)
  ## then the inequalities
  if (Nineq > 0) {
    con <- rbind(con, G)
    rhs <- c(rhs, H)
    dir <- c(dir, rep(">=", Nineq))
  }
  
  AllX   <- NULL
  Summed <- rep(0,Nx)
  nsum   <- 0

  if (ispos) {

    for (i in 1:Nx) {
      obj <- rep(0, Nx)
      obj[i] <- 1
      lmin <- lp("min", obj, con, dir, rhs)
      if(lmin$status == 0) Range[i, 1] <- lmin$objval else
      if(lmin$status == 3) Range[i, 1] <- -1e30 else
      Range[i, 1] <- NA
      lmax <- lp("max", obj, con, dir, rhs)
      if(lmax$status == 0) Range[i, 2] <- lmax$objval else
      if(lmax$status == 3) Range[i, 2] <- 1e30 else
      Range[i, 2] <- NA
      if (central) {
        if (! any (is.na(lmin$solution)) && lmin$status==0) {
          Summed<- Summed + lmin$solution
          nsum<-nsum+1
        }
        if (! any (is.na(lmax$solution)) && lmax$status==0) {
          Summed<- Summed + lmax$solution
          nsum<-nsum+1
        }
      }
      if (full) {
        if( lmin$status==0)     AllX <- cbind(AllX,lmin$solution)
        if( lmax$status==0)     AllX <- cbind(AllX,lmax$solution)
      }
    }
  } else{
    ## First test if problem is solvable...
    Sol <- lsei(E=E,F=F,G=G,H=H)
    if (Sol$residualNorm > tol)  {
      warning (paste("cannot proceed: problem not solvable at requested tolerance",tol))
      return(Range)
    }

    ## double the number of unknowns: x -> x1 -x2, where x1>0 and x2>0
    con <- cbind(con,-1*con)

    for (i in 1:Nx) {
      obj <- rep(0, 2*Nx)
      obj[i]    <- 1
      obj[Nx+i] <- -1

      lmin <- lp("min", obj, con, dir, rhs)
      if(lmin$status == 0) Range[i, 1] <- lmin$objval else
      if(lmin$status == 3) Range[i, 1] <- -1e30 else
      Range[i, 1] <- NA
      lmax <- lp("max", obj, con, dir, rhs)
      if(lmax$status == 0) Range[i, 2] <- lmax$objval else
      if(lmax$status == 3) Range[i, 2] <- 1e30 else
      Range[i, 2] <- NA
      if (central) {
        if (! any (is.na(lmin$solution)) && lmin$status==0) {
          Summed<- Summed + lmin$solution
          nsum<-nsum+1
        }
        if (! any (is.na(lmax$solution)) && lmax$status==0) {
          Summed<- Summed + lmax$solution
          nsum<-nsum+1
        }
      }
      if (full) {
        if( lmin$status==0)     AllX <- cbind(AllX,lmin$solution)
        if( lmax$status==0)     AllX <- cbind(AllX,lmax$solution)
      }
    }
    if (central) Summed <- Summed[1:Nx]- Summed[(Nx+1):(2*Nx)]
    if (full)    AllX   <- AllX[1:Nx,]-AllX[(Nx+1):(2*Nx),]
  }
  colnames(Range) <- c("min", "max")
  xnames <- colnames(E)
  if (is.null(xnames))
    xnames <- colnames(G)
  rownames(Range) <- xnames
  if (central) Range<-cbind(Range, central = Summed/nsum)
  if (full) Range<-cbind(Range, AllX)
  return(Range)
}
