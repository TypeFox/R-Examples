Hreprrational <- function(pts) { ## input numeric or rational, output rational
  #print(nrow)
  colNames <- colnames(pts)
  if(class(pts[1])=="character") {
    foo <- cbind("0", "1", pts)
  } else {
    # it's important to use rationals in scdd.
    ###########   nr <- nrow(pts)
    foo <- cbind(0, 1, pts)
    foo <- d2q(foo)
  }
  foorational <- scdd(foo, representation = "V", roworder="maxcutoff")$output # H repr from V repr.
  if (!is.null(colNames)) colnames(foorational) <- c("eq", "b", colNames) ## not very useful ?
  return(foorational)
}
