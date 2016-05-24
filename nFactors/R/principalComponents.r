"principalComponents" <-
function(R) {
 nVar   <- dim(R)[2]
 acp             <- eigen(R)
 values          <- acp$values
 vectors         <- acp$vectors # Normed vecteurs to 1
 varExplained    <- round((values/nVar)*100,    2)
 cumVarExplained <- round(cumsum(varExplained), 2)
 loadings        <- vectors  %*% diag(values^0.5)  # F1 * diag(E)
 acp             <- list(values          = values,
                         varExplained    = varExplained,
                         cumVarExplained = cumVarExplained,
                         vectors         = vectors,
                         loadings        = loadings)
 return(acp)
 }

