"mhalf" <-
function(A) {
   d <- diag(eigen(A)$values)
   v <- eigen(A)$vectors
   Ah <- v%*%solve(sqrt(d))%*%t(v)
}

