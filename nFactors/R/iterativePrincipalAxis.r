iterativePrincipalAxis <-
function(R, nFactors=2, communalities="component", iterations=20, tolerance=0.001) {
 if (communalities == "component")            diag(R)  <- componentAxis(R)$communalities
 if (communalities == "maxr")      { RT <- R; diag(RT) <- 0; diag(R) <- apply(RT, 1, max)}
 if (communalities == "ginv")                 diag(R)  <- sqrt(1-1/diag(ginv(R)))
 if (communalities == "multiple")  {
  if (all(eigen(R)$values > 0)) diag(R) <- sqrt(1-1/diag(solve(R)))  # Gorsuch (1983, p. 106)
  else return("Not all eigenvalues are grater than 0") # Verication of positive definiteness
  }
  iter <- 1; tol <- 1
  while ((iter < iterations) && (tol > tolerance)) {     # for (i in (1:iterations))
   oldR    <- diag(R)
   diag(R) <- componentAxis(R, nFactors)$communalities
   tol     <- max(abs(diag(R) - oldR))
   iter    <- iter + 1
  }
 if (tol > tolerance) warning("Maximum number of iterations needed before the desired tolerance: cautious solution.")
 iapa <- componentAxis(R, nFactors)
 iapa <- list(values          = iapa$values,
              varExplained    = iapa$varExplained,
              cumVarExplained = iapa$cumVarExplained,
              loadings        = iapa$loadings,
              iterations      = iter,
              tolerance       = tol)
 return(iapa)
 }

