psdcor <- function(m, method = c("higham", "eigen"), ...) {
  method <- match.arg(method)
  m <- as.matrix(m)
  if (length(m) <= 1) {
    ret <- m
  } else {
    ret <- switch(method,
      eigen = {
        require(sfsmisc)
        posdefify(m, ...)
      },
      higham = {
        require(sfsmisc)
        nearcor(m, ...)$cor
  #      require(Matrix);
  #      as.matrix(nearPD(m, corr=TRUE, ...)$mat)
      }
    )
  }
  ret
}
