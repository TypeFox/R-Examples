#Jacobian for ordered model
jacobian <- function(kappa){
  ## ij elemenet is dkappa_i/d alpha_j  
  k <- length(kappa)
  etheta <- exp(kappa)
  mat <- matrix(0 ,k, k)
  for (i in 1:k) mat[i:k, i] <- etheta[i]
  mat
}

# Compute time
make.time <- function(object){
  et <- object$time[3]
  s <- round(et,0)
  h <- s%/%3600
  s <- s-3600*h
  m <- s%/%60
  s <- s-60*m
  paste(h, "h:", m, "m:", s, "s", sep="")
}

#' @export
ordinal <- function(link = c('probit', 'logit')){
  link <- match.arg(link)
  list(family = 'ordinal', link = link)
}

# Added in version 0.2
repRows <- function(x, n){
  matrix(rep(x, each = n), nrow = n)
}

# Added in version 0.2
make.add <- function(row, col, Ka){
  sa <- makeL(1:rep(0.5 * Ka * (Ka + 1)))
  for (k in row:col){ 
    cb <- sa[k, row]
    form <- paste(paste("x",  cb:(cb + (Ka - k)), sep = ""), paste("x", cb, sep = ""), sep = "*")
  }
  form
}
