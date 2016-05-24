npmr <-
function(X, Y, lambda = exp(seq(7, -2)), s = 0.1, eps = 1e-6,
    group = NULL, accelerated = FALSE, B.init = NULL, b.init = NULL,
    quiet = TRUE) {

    if (is.null(dim(Y))) {
        colnames = sort(unique(Y))
        Y = model.matrix(~ as.factor(Y) - 1)
        colnames(Y) = colnames
    }

    if (is.null(B.init)) {
#        B = matrix(rnorm(ncol(X)*ncol(Y)), ncol(X), ncol(Y))
        B = matrix(0, ncol(X), ncol(Y))
    } else B = B.init

    if (is.null(b.init)) {
        b = log(colMeans(Y)) - mean(log(colMeans(Y)))
    } else b = b.init

    if (nrow(X) != nrow(Y)) {
        stop('X and Y do not have matching numbers of observations')
    } else if (!is.null(group) & length(group) != ncol(X)) {
        stop('X and group do not have matching numbers of observations')
    } else if (ncol(X) != nrow(B)) {
        stop('Number of rows of B.init does not match number of columns of X')
    } else if (ncol(Y) != ncol(B)) {
        stop('Number of columns of B.init does not match number of classes')
    } else if (ncol(Y) != length(b)) {
        stop('Length of b.init does not match number of classes')
    }

    B.path = array(NA, dim = c(ncol(X), ncol(Y), length(lambda)))
    b.path = array(NA, dim = c(ncol(Y), length(lambda)))
    objective.path = rep(NA, length(lambda))

    cat('Progress: ')
    for (l in 1:length(lambda)) {
        solution = PGDnpmr(B, b, X, Y, lambda[l], s = s, group = group,
            accelerated = accelerated, eps = eps, quiet = quiet)
        B = solution$B
        b = solution$b
        B.path[, , l] = as.matrix(B)
        b.path[, l] = scale(b, scale = FALSE)
        objective.path[l] = min(solution$objectivePath)
        cat(round(100*l/length(lambda)))
        cat('% ')
    }
    cat('\n')

    rownames(B.path) = colnames(X)
    colnames(B.path) = rownames(b.path) = colnames(Y)

    fit = list(call = match.call(), B = B.path, b = b.path,
        objective = objective.path, lambda = lambda)
    class(fit) = 'npmr'
    return(fit)
}
