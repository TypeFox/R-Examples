smirnov <-
function(x) {
x <- as.matrix(x)
x <- (x > 0) * 1
x <- x[, !(colSums(x) == 0 | colSums(x) == nrow(x))]
s <- nrow(x)
n <- ncol(x)
w.0 <- colSums(x == 1)/colSums(x == 0)
w.1 <- 1/w.0
mat <- matrix(0, s, s)
for (i in 1:s) {
        for (j in 1:s) {
        ini <- x[i, ] + x[j, ]
        txy.0 <- (ini == 0) * w.0
        txy.1 <- (ini == 2) * w.1
        txy.r <- (ini == 1)
        txy <- 1/n * sum(txy.0 + txy.1 - txy.r)
        mat[i, j] <- txy
        }
}
mat.scaled <- 2*(mat+1)/s
diag(mat.scaled) <- (diag(mat)*(s-1)-1)/(s*(s-2))
mat.scaled
}

