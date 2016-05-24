pathCoeff <- function(object, ...){
  UseMethod("pathCoeff", object)
}


# Calculates the matrix of path coefficients.
pathCoeff.default <- function(model, factor_scores, method, pairwise, ...){
  ifelse(pairwise, use <- "pairwise.complete.obs", use <- "everything")
  latent <- model$latent             # names of latent variables
  strucmod <- model$strucmod         # names of manifest variables
  fscores <- factor_scores           # matrix of estimated factor scores
  pC <- matrix(0, nrow=length(latent), ncol=length(latent))
  rownames(pC) <- latent
  colnames(pC) <- latent

  for(i in latent){
    if (i %in% strucmod[,2]){
      # which latents are direct predecessors from latent[i]
      # latent[i] is the dependent variable
      index <- which(i == strucmod[,2])
      indpnt <- strucmod[index,1] # the independent LVs

      # solving the structural equation for latent[i]
      #pC[indpnt, i] <- solve(cor(as.matrix(fscores[,indpnt]), use=use, method=method)) %*%
      #                 cor(fscores[,indpnt], fscores[,i], use=use, method=method)
      pC[indpnt, i] <- solve(cor(fscores[,indpnt, drop=FALSE], use=use, method=method),
                       cor(fscores[,indpnt], fscores[,i], use=use, method=method))
    }
  }
  return(pC)
}

pathCoeff.sempls <- function(object, ...){
  coeffs <- object$path_coefficients
  class(coeffs) <- "pathCoeff"
  return(coeffs)
}

print.pathCoeff <- function(x, na.print=".", digits=2, abbreviate=FALSE, ...){
  pathCoeff <- x
  pathCoeff[pathCoeff==0] <- NA
  if(abbreviate) dimnames(pathCoeff) <- lapply(dimnames(pathCoeff), abbreviate, ...)
  print.table(pathCoeff, na.print=na.print, digits=digits, ...)
  invisible(x)
}
  
