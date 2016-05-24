##==============================================================================
## lsei= least squares with equality and inequality constraints
##==============================================================================

Lsei <- function(...) UseMethod ("Lsei")
Lsei.character <- function(...) Lsei.limfile(...)
Lsei.double <- function(...) lsei(...)

##==============================================================================
# reads an input file and solves the linear inverse model using lsei
##==============================================================================
Lsei.limfile <- function(file, exact =NULL, parsimonious = FALSE,
         verbose=TRUE, ...) {

  # read file and create inverse matrices
  lim    <- Setup.limfile(file,verbose=verbose)

  Lsei.lim(lim,exact,parsimonious=parsimonious,...)      # solve
}

##==============================================================================
## Solves inverse model, using lsei
## lim is a list that contains the input matrices
## if parsimonious is true: sum of squared flows has to be minimised too
##==============================================================================

Lsei.lim <- function(lim, exact =NULL, parsimonious = FALSE, ...) {
# 0. Setup problem

  Ncomp  <- lim$NComponents
  Nx     <- lim$NUnknowns

  A <- E <- NULL
  B <- F <- NULL
  Napp <- nrow(lim$A)

  if (is.null(exact)) {   # Equalities: all equations are approximate
    Neq  <- 0
    E      <- lim$A[]
    F      <- lim$B[]

  } else  {               # Equalities and approximate equations
    if (max(exact) > Napp)
      stop("error: cannot solve Lsei.lim. Equations to be met exactly do not exist")
    Neq  <- length(exact)
    E    <- lim$A[exact,]
    F    <- lim$B[exact]
    # Approximate equations : MIN |AX-B|

    Napp   <- Napp-Neq    # Number approximate
    if (Napp > 0) {
      A      <- lim$A[-exact,]
      B      <- lim$B[-exact]
    }
  }

 if (parsimonious)   {     # min sum of squared flows
   A <- rbind(A,diag(nrow=Nx,ncol=Nx))
   B <- c(B,rep(0,Nx))
 }

 if (is.null (A) )
   warning("Lsei.lim - there are no approximate equations in lsei!")

  # Inequalities G*X>=H
  G <- lim$G
  H <- lim$H


  sol <- lsei(A,B,E,F,G,H,...)
  names(sol$X) <- lim$Unknowns
  if (sol$IsError)
    warning("Problem could not be solved - least squares solution returned")

  return(sol)
}
