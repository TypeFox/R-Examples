print.npmr <-
function(x, ...) {

    rank = rep(NA, length(x$lambda))

    for (l in 1:length(x$lambda)) {
        rank[l] = qr(x$B[,,l])$rank
    }

    cat('\n Call:\n')
    print(x$call)
    cat('\n')
    print(data.frame(lambda = x$lambda, rank = rank, objective = x$objective))
    cat('\n')
}
