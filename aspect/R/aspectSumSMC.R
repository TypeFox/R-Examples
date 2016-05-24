`aspectSumSMC` <-
function(r) {
  m <- dim(r)[1]
  f <- 0
  g <- matrix(0,m,m)
  for (j in 1:m) {
    a <- aspectSMC(r,j)
    f <- f+(a$f)
    g <- g+(a$g)
  }
  list(f = f, g = g)
}

