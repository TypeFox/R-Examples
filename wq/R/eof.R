eof <-
function (x, n, scale. = TRUE) {

  # Validate args
  if (!is.matrix(x) && !is.data.frame(x))
    stop("x must be a 'matrix' or 'data.frame'")
  if (identical(colnames(x), NULL))
    colnames(x) <- paste("v", 1:ncol(x), sep="")
  if (anyDuplicated(colnames(x)) > 0)
    stop("x must have distinct column names")
  if (is.mts(x))
    rownames(x) <- time(x)
  if (identical(rownames(x), NULL))
    rownames(x) <- 1:nrow(x)
  if (anyDuplicated(rownames(x)) > 0)
    stop("x must have distinct row names")

  # get EOFs (as scaled eigenvectors)
  pr1 <- prcomp(x, scale.=TRUE)
  eigenval1 <- pr1[["sdev"]][1:n]^2
  eigenvec1 <- pr1[["rotation"]][, 1:n]
  eof1 <- eigenvec1 %*% diag(sqrt(eigenval1), n, n)
  scores1 <- pr1[["x"]][, 1:n]
  amp1 <- scale(scores1)

  # get REOFs by orthogonally rotating EOFs
  if (identical(n, 1)) {
    reof <- as.matrix(eof1)
    amp <- amp1
  } else {
    pr2 <- varimax(eof1)
    reof <- unclass(pr2[["loadings"]])
    rotater <- pr2[["rotmat"]]
    amp <- amp1 %*% rotater
  }

  # Eigenvalues and total percent variance for first n
  eigs <- pr1[["sdev"]]^2
  eigen.pct <- round(100 * eigs/sum(eigs), 1)
  totvar.pct <- round(100 * cumsum(eigs/sum(eigs)), 1)

  # return results
  colnames(reof) <- colnames(amp) <- paste('EOF', 1:n, sep='')
  reof <- cbind(id=ordered(colnames(x), levels=colnames(x)),
  	as.data.frame(reof))
  amp <- cbind(id=ordered(rownames(x), levels=rownames(x)),
  	as.data.frame(amp))
  rownames(reof) <- rownames(amp) <- NULL
  list(REOF=reof, amplitude=amp, eigen.pct=eigen.pct,
  	variance=totvar.pct)
}
