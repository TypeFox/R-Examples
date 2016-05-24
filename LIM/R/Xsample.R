##==============================================================================
##  Samples a linear inverse problem
##==============================================================================

Xsample <- function(lim, exact =NULL, ...) {
  ## 0. Setup problem

  Ncomp  <- lim$NComponents
  Nx     <- lim$NUnknowns

  A <- E <- NULL
  B <- F <- NULL
  Napp <- nrow(lim$A)

  if ( is.null(exact) )  { # Equalities: all equations are approximate
    Neq  <- 0
    E    <- lim$A[]
    F    <- lim$B[]

  } else  {               # Equalities and approximate equations
    if ( max(exact) > Napp )
      stop("error: cannot solve Lsei.lim. Equations to be met exactly do not exist")
    Neq  <- length(exact)
    E    <- lim$A[exact,]
    F    <- lim$B[exact]
    # Approximate equations : MIN |AX-B|

    Napp   <- Napp-Neq    # Number approximate

    if (Napp > 0)  {
      A      <- lim$A[-exact,]
      B      <- lim$B[-exact]
    }
  }

  # Inequalities G*X>=H
  G <- lim$G
  H <- lim$H

  res   <- xsample (A=A, B=B, E=E, F=F,
                    G=G, H=H, ...)$X
  colnames(res) <- lim$Unknowns

  return(res)

}
