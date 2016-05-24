"half" <-
function(A) {
   d <- diag(eigen(A)$values)
   v <- eigen(A)$vectors
   Ah <- v%*%sqrt(d)%*%t(v)
}

