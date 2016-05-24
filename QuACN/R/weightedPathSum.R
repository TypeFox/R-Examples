.weightedPathSum <- function(e, path, weightfunc, unit=1, rmdupls=TRUE) {
  if (is.null(path))
    sum(sapply(names(e), function(v)
      .weightedPathSum(e, v, weightfunc, unit=unit, rmdupls=rmdupls)))
  else {
    corr <- if (rmdupls) 0.5 else 1

    if (length(path) == 1)
      result <- unit
    else
      result <- prod(sapply(1:(length(path) - 1), function(i) {
        weightfunc(i, path[[i]], path[[i + 1]])
      })) * corr

    last <- path[[length(path)]]
    result + sum(sapply(e[[last]], function(v) {
      if (v %in% path)
        0
      else
        .weightedPathSum(e, c(path, v), weightfunc, unit=unit, rmdupls=rmdupls)
    }))
  }
}
