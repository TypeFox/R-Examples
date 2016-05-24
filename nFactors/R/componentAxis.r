"componentAxis" <-
function(R, nFactors=2) {
  nVar            <- dim(R)[2]
  acp             <- principalComponents(R)
  values          <- acp$values[(1:nFactors)]
  varExplained    <- round((values/nVar)*100,    2)
  cumVarExplained <- round(cumsum(varExplained), 2)
  loadings        <- acp$vectors[,(1:nFactors)]  %*% diag(values^0.5)  # F1 * diag(E)
  communalities   <- apply(loadings*loadings,1,sum)
  apa             <- list(values          = values,
                          varExplained    = varExplained,
                          cumVarExplained = cumVarExplained,
                          loadings        = loadings,
                          communalities   = communalities)
  return(apa)
  }
