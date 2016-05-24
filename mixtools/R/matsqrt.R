matsqrt <- function(x) {
xe <- eigen(x)
xe.v <- xe$values
if(all(xe.v >= 0)) {
xe.v1 <- diag(sqrt(xe.v))
}
xvalue <- cbind(xe$vectors)
xvalue.1 <- solve(xvalue)
out <- xvalue %*% xe.v1 %*% xvalue.1
out
}
