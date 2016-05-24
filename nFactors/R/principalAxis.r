"principalAxis" <-
function(R, nFactors=2, communalities="component") {
 if (communalities == "component")            diag(R)  <- componentAxis(R)$communalities
 if (communalities == "maxr")      { RT <- R; diag(RT) <- 0; diag(R) <- apply(RT, 1, max)}
 if (communalities == "ginv")                 diag(R)  <- sqrt(1-1/diag(ginv(R)))
 if (communalities == "multiple")  {
  if (all(eigen(R)$values > 0)) diag(R) <- sqrt(1-1/diag(solve(R)))  # Gorsuch (1983, p. 106)
  else return("Not all eigenvalues are greater than 0") # Verication of positive definiteness
  }
 apa <- componentAxis(R, nFactors)
 return(apa)
 }

