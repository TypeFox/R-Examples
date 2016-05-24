print.cv.npmr <-
function(x, ...) {

    rank = qr(x$fit$B[,,1])$rank
    error = x$error[which(x$lambda == x$lambda.min)]/x$n

    cat('\n Call:\n')
    print(x$call)
    cat('\n')
    cat('Optimal lambda: ')
    cat(x$lambda.min, '\n')
    cat('Rank of fitted coefficient matrix: ')
    cat(rank, '\n')
    cat('Mean CV error: ')
    cat(error, '\n\n')
}
