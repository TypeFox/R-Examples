`bmat` <-
function(diss, wgths, d, eps = 1E-12)
{
  z <- ifelse(d < eps, 1, 0)
  b <- as.matrix((wgths*diss*(1-z))/(d+z))
  r <- rowSums(b) 
  return(diag(r)-b)
}

