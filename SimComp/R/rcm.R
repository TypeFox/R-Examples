rcm <-
function(nrow=NULL, ncol=NULL) {


if (!is.null(nrow) & !is.null(ncol)) {
  stop("Either nrow or ncol must be specified")
}
if (is.null(nrow)) {
  if (is.null(ncol)) {
    stop("Either nrow or ncol must be specified")
  } else {
    nrow <- ncol
  }
}

if (!is.numeric(nrow) | length(nrow)!=1) {
  stop("nrow or ncol must be a single numeric value")
}

mat <- rmvnorm(nrow,rep(0,nrow),diag(nrow))
for (i in 1:nrow) {
  mat[i,] <- mat[i,]/sqrt(t(mat[i,])%*%mat[i,])
}
return(mat%*%t(mat))


}
