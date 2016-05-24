predict.npmr <-
function(object, newx, ...) {

    if (ncol(newx) != nrow(object$B)) {
        stop('Number of variables in newx does not match model fit')
    }

    nlambda = ncol(object$b)
    eta = P = array(NA, c(nrow(newx), dim(object$B)[2], nlambda))
    for (l in 1:nlambda) {
        eta[, , l] = as.matrix(matrix(1, nrow(newx), 1) %*% t(object$b[, l]) +
            newx %*% object$B[, , l])
        P[, , l] = exp(eta[, , l])/rowSums(exp(eta[, , l]))
    }
    P
}
