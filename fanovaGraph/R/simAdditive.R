simAdditive <- function(newdata, mu, parameter, covtype, cl, iso=FALSE, eps.R = 1e-08){
  newdata <- as.matrix(newdata)
  n <- nrow(newdata)
  R <- matrix(0, ncol = n, nrow = n)
  ncl <- length(cl)
  if (length(iso)==1) 
    iso <- rep(iso,ncl)
  parameter <- paramList2Vect(parameter, cl, iso)
  alpha <- parameter[1:ncl]
  theta <- parameter[-(1:ncl)]
  thetalist <- cl
  # building thetalist:
  ntemp <- 0
  nclani <- ncl - sum(iso)
  # for (j in (1:nclani)) { # DOES NOT RESPECT ORDER
  for (j in (1:ncl)[!iso]){
    thetalist[[j]] <- theta[(ntemp + 1):(ntemp + length(cl[[j]]))]
    ntemp <- ntemp + length(cl[[j]])
  }
  if (nclani < ncl) {
  # for (j in ((nclani + 1):ncl)) { # DOES NOT RESPECT ORDER
    for (j in (1:ncl)[iso]) {
      thetalist[[j]] <- theta[ntemp]
      ntemp <- ntemp + 1
    }
  }
  #building R
  for (j in 1:ncl) { 
    cor.str <- covStruct.create(covtype = covtype, d = length(cl[[j]]), 
            var.names = NULL, known.covparam = "All", coef.cov = thetalist[[j]], 
            coef.var = alpha[j], iso = iso[j])
    R <- R + covMatrix(object = cor.str, X = as.matrix(newdata[, cl[[j]]]))[[1]]
  }
  R <- R + diag(eps.R, n, n)
  whitenoise <- rnorm(n)
  cR <- chol(R)
  y <- mu + t(cR) %*% whitenoise   # t() because unlike Rasmussen cR'cR=R
  return(drop(y))
}
