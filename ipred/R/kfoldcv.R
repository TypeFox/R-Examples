# $Id: kfoldcv.R,v 1.3 2002/09/12 08:56:42 hothorn Exp $

kfoldcv <- function(k,N, nlevel=NULL) {
  if (is.null(nlevel)) {
    # no stratification
    if (k > N) return(c(rep(1, N), rep(0, k-N)))
    fl <- floor(N/k)
    ce <- ceiling(N/k)
    if (fl == ce) return(rep(fl, k)) 
      else 
    return(c(rep(ce, round((N/k - fl)*k)), rep(fl, round((1 - (N/k -
                     fl))*k))))
  } else {
    # stratification
    # if (!is.integer(nlevel)) stop("nlevel is not a vector if integers")
    kmat <- matrix(0, ncol=k, nrow=length(nlevel))
    for (i in 1:length(nlevel))
      kmat[i,] <- kfoldcv(k, nlevel[i])
    return(kmat)
  }
}
