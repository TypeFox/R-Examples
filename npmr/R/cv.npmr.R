cv.npmr <-
function(X, Y, lambda = exp(seq(7, -2)), s = 0.1, eps = 1e-6,
    group = NULL, accelerated = FALSE, B.init = NULL, b.init = NULL,
    foldid = NULL, nfolds = 10) {

    if (is.null(dim(Y))) {
        colnames = sort(unique(Y))
        Y = model.matrix(~ as.factor(Y) - 1)
        colnames(Y) = colnames
    }

    if (is.null(foldid)) foldid = sample(rep(1:nfolds, length = nrow(X)))
    nfolds = length(unique(foldid))

    if (nrow(X) != nrow(Y)) {
        stop('X and Y do not have matching numbers of observations')
    } else if (nrow(X) != length(foldid)) {
        stop('X and foldid do not have matching numbers of observations')
    }

    error = matrix(NA, nfolds, length(lambda))

    for (fold in unique(foldid)) {
        fit = npmr(X[foldid != fold, ], Y[foldid != fold, ], lambda,
            s = s, eps = eps, group = group, accelerated = accelerated,
            B.init = B.init, b.init = b.init)
        for (l in 1:length(lambda)) {
            error[sort(unique(foldid)) == fold, l] =
                sum(-log(rowSums(predict(fit, X[foldid == fold, ])[, , l]*
                Y[foldid == fold, ])))
        }
    }

    lambda.min = lambda[which.min(colSums(error))]
    if (lambda.min == min(lambda)) {
        warning('lambda chosen through CV is smallest lambda')
    } else if (lambda.min == max(lambda)) {
        warning('lambda chosen through CV is largest lambda')
    }
        
    fit = list(call = match.call(), error = colSums(error),
        fit = npmr(X, Y, lambda.min), lambda.min = lambda.min, lambda = lambda,
        n = nrow(X))
    class(fit) = "cv.npmr"
    return(fit)
}
