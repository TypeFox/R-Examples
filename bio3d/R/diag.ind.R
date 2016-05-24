`diag.ind` <-
function (x, n=1, diag = TRUE) {
  x <- as.matrix(x)
  if (diag) {
    !(row(x) > col(x)) + (row(x) <= col(x)-n)
  } else { !(row(x) >= col(x)) + (row(x) <= col(x)-n) }
}

