colVarsC <- function(x)
{
  if (!is.matrix(x)) x <- as.matrix(x)
  CN <- colnames(x)
  if (NCOL(x) == 1) stop("'x' should have more than one column!")
  OUT <- .Call("colVarsC", x, PACKAGE = "propagate")
  names(OUT) <- CN
  OUT
}

rowVarsC <- function(x)
{
  if (!is.matrix(x)) x <- as.matrix(x)
  RN <- rownames(x)
  if (NROW(t(x)) == 1) stop("'x' should have more than one row!")
  OUT <- .Call("rowVarsC", x, PACKAGE = "propagate")
  names(OUT) <- RN
  OUT
}
