miqr.solve <- function(a, b) {
  x <- tryCatch({
    qr.solve(a, b);
  }, error = function(e) {
    e;
  }, finally = "");

  if (!is.matrix(x)) {
    x <- ginv(a) %*% b;
  }

  x;
} # miqr.solve()
