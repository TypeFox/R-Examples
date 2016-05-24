logL <-
function(B, b, X, Y) {
    eta = matrix(1, nrow(X), 1) %*% t(b) + X %*% B
    sum(Y*(eta - log(rowSums(exp(eta)))))
}
