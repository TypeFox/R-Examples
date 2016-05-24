##==============================================================================
## Linear programming
##==============================================================================

Linp <- function(...) UseMethod ("Linp")
Linp.character <- function(...) Linp.limfile(...)
Linp.double <- function(...) linp(...)

##==============================================================================
# reads an input file and solves the linear programming inverse model
##==============================================================================

Linp.limfile <- function(file,verbose=TRUE,...)  {
  lim    <- Setup.limfile(file,verbose=verbose)
  Linp.lim(lim,...)
}

##==============================================================================
## Linp.lim : Solves linear programming model
##==============================================================================

Linp.lim <- function(lim,cost=NULL,ispos=lim$ispos,...) {

  Nx <- lim$NUnknowns
  B  <-lim$B
  H  <-lim$H
  A  <- NULL
  G  <- NULL

  if ( is.null(cost))
    Cost   <- lim$Cost
  else Cost <- cost

  if ( is.null(cost))
    Profit <- lim$Profit
  else Profit <- NULL
  
  if (! is.null(Cost  )&& is.vector(Cost  ))
    Cost   <- matrix(data=Cost, nrow=1)
  if (! is.null(Profit)&& is.vector(Profit))
    Profit <- matrix(data=Profit, nrow=1)

  if (ispos) {
    A <- lim$A
    G <-lim$G
  } else {
    if (! is.null(lim$A))
      A <- cbind(lim$A,-1*lim$A)
    if (! is.null(lim$G))
      G <- cbind(lim$G, -1*lim$G)
    if (! is.null(Cost))
      Cost   <- cbind(Cost, -1 * Cost)
    if (! is.null(Profit))
      Profit <- cbind(Profit, -1 * Profit)
  }
  res <- NULL
  i1 <- 0
  if (!is.null(Cost)) {
    for ( i in 1:nrow(Cost)) {
      resi <- linp (A,B,G,H,Cost[i,],...)
      res$residualNorm <- c(res$residualNorm, resi$residualNorm)
      res$solutionNorm <- c(res$solutionNorm, resi$solutionNorm)
      res$X <- rbind(res$X,matrix(data=resi$X, nrow=1))
      rownames(res$X)[i]<-as.character(lim$costnames[i])
    }
    i1 <- i
  } else if (is.null(lim$Profit)) return(NULL)

  if (! is.null(Profit)) {
    for ( i in 1:nrow(Profit))  {
      resi<-linp (A,B,G,H,-1*Profit[i,],...)
      res$residualNorm <- c(res$residualNorm, resi$residualNorm)
      res$solutionNorm <- c(res$solutionNorm, -resi$solutionNorm)
      res$X <- rbind(res$X,matrix(data=resi$X,nrow=1))
      rownames(res$X)[i1+i]<-as.character(lim$profitnames[i])
    }
  }
  if (!ispos)
    res$X <-  res$X[,1:Nx]-res$X[,(Nx+1):(2*Nx)]
  res$X <- matrix(data=res$X, ncol=Nx)
  colnames(res$X) <- lim$Unknowns

  return(res)
}
