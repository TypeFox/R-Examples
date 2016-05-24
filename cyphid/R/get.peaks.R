get.peaks <-
function (x, span = 40) 
  {
    z <- embed(rev(as.vector(x)), dimension = span)
    z <- z[rev(seq(nrow(z))), ]
    s <- span%/%2
    v <- max.col(z, ties.method = "first") == 1 + s
    z <- c(rep(FALSE, s), v)
    ans <- c(z[1:(length(z) - s)], rep(FALSE, span - 1))
    nx <- NCOL(x)
    if (nx > 1) 
      matrix(ans, ncol = nx)
    else ans
  }

