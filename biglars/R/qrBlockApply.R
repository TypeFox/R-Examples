
qrBlockApply <- 
function(x, y = NULL, blockSize = NULL)
{
  ##
  ## block QR factorization of matrix x
  ##
  nrowx <- nrow(x)
  ncolx <- ncol(x)

  if(!is.ff(x))
    x <- as.ff(x)
  
  if(is.null(y)) {
    y <- numeric(nrowx)
  }
  if(!is.ff(y)) {
    y <- as.ff(y)
  }


  ## "ffbatchbytes" curently defaults to 2^16.
  #
  if(is.null(blockSize)) {
    blockSize <- getOption("ffbatchbytes")
  }
  else if(blockSize < ncolx) {
    blockSize <- ncolx
    warning(paste("block size has been increased to", ncolx, sep = " "))
  }

  # 'i1', 'i2' are the bounds of the successive chunks, determined by
  # ffrowapply().
  #
  # 'fit' is a structure consisting of a Cholesky factor, 'R', and its
  # effects vector, 'Qty'.
  #
  i1 <- i2 <- 0L
  dimR <- ncolx + 1
  fit <- list(R = matrix(0, dimR, dimR), Qty = as.double(rep(0, dimR)))
  
  ffrowapply( {
    rowRange <- i1:i2
    fit <- qrBlock(x[rowRange, , drop=FALSE], y[rowRange], fit$R, fit$Qty)
    }, X=x, BATCHSIZE = blockSize)

  R <- fit$R
  Qty <- fit$Qty
  R[row(R) > col(R)] <- 0

  namX <- dimnames(x)[[2]]
  if(is.null(namX))
    namX <- 1:ncolx
  dimnames(R) <- list(NULL, c("(Intercept)", namX))
  names(Qty) <- NULL

  out <- list(R, Qty)
  names(out) <- c("R", "Qty")
  out
}

