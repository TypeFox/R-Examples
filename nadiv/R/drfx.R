drfx <- function(G, fac, dataf, ...){
   dataf[, fac] <- as.factor(dataf[, fac])
   d <- nrow(G)
   if(all(G == G[1,1]) & d > 1){
      warning("variance-covariance matrix 'G' may have caused 'chol.default(G)' error.  If so, consider subtracting 0.0001 from the covariances to make correlations < 1 or >-1")
   }
   Z <- sparse.model.matrix(as.formula(paste0("~", fac, " - 1")), dataf)
   M <- grfx(n = ncol(Z), G = G, incidence = Diagonal(ncol(Z)), saveIncidence = FALSE, ...)
   fx <- sapply(seq.int(d), FUN = function(c){ (Z %*% M[, c])@x})
 return(list(fx = fx, Z = Z))
}

