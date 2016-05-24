

print.rvfactor <- function(x, all.levels=FALSE, ...) {
  s <- summary(x, all.levels=all.levels)
  ds <- dimnames(s)
  if (!is.null(.names <- names(x))) {
    s <- cbind(name=.names, s)
  } else if (!is.null(.dn <- dimnames(x))) {
    sud <- rvpar("summary.dimnames")
    if (!is.null(sud) && !is.na(sud) && is.logical(sud)) {
      da <- lapply(.dn, function (na) if (is.null(na)) rep("", nrow(s)) else na)
      rw <- da[[1]][row(x)]
      cl <- da[[2]][col(x)]
      s <- cbind(row=rw, col=cl, " "=":", s)
    }
  }
  print(s)
}

