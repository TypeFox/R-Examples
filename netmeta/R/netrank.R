netrank <- function(x, small.values="good"){
  
  ## Check for netmeta object
  ##
  meta:::chkclass(x, "netmeta")
  small.values <- meta:::setchar(small.values, c("good", "bad"))
  
  
  ## Calculate one-sided p-values
  ##
  w <- (1 + sign(x$TE.random))/2
  p <- x$pval.random
  ##
  if (small.values=="good")
    P <- w*p/2 + (1-w)*(1 - p/2)
  else
    P <- w*(1 - p/2) + (1-w)*p/2
  
  
  ## Row means provide P-scores
  ##
  Pscore <- rowMeans(P[,], na.rm=TRUE)
  
  
  res <- list(Pscore=Pscore,
              Pmatrix=P,
              small.values=small.values,
              x=x,
              title = x$title,
              version=packageDescription("netmeta")$Version)

  class(res) <- "netrank"
  
  res
}
