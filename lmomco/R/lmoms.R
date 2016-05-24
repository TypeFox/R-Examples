"lmoms" <-
function(x,nmom=5) {
  n <- length(x)

  if(nmom > n) {
    stop("More L-moments requested by parameter 'nmom' than data points available in 'x'")
  }

  if(length(unique(x)) == 1) stop("all values are equal--Lmoments can not be computed")
  z <- TLmoms(x,nmom=nmom)
  z$source <- "lmoms"
  return(z)
}
