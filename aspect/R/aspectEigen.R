`aspectEigen` <-
function(r, p = 1) {
  s <- eigen(r)
  v <- s$vectors[,1:p]
  list(f = sum(s$values[1:p]), g = v%*%t(v))
}

