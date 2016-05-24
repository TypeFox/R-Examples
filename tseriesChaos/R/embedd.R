embedd <- function (x, m, d, lags) {
  if (is.null(colnames(x))) {
    vars <- sprintf("V%d", 1:NCOL(x))
  } else {
    vars <- colnames(x)
  }

  x <- as.ts(x)
  if (missing(lags)) {
    checkEmbParms(x, m, d)
    lags <- ((1:m) - 1) * d
  }

  names <- sprintf("%s/%d", vars, lags[1])
  res <- lag(x, lags[1])
  if(length(lags) > 1) {
    for (i in 2:length(lags)) {
      res <- ts.intersect(res, lag(x, lags[i]))
      names <- append(names, sprintf("%s/%d", vars, lags[i]))
    }
  }

  res <- matrix(res, nrow = nrow(res), ncol = ncol(res))
  colnames(res) <- names
  res
}
