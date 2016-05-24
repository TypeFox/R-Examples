if (getRversion() >= "2.15.1") {
  # These are available in the suggested packages:
  # penalized & nnet
  utils::globalVariables(c("penalized", "multinom"));
}

slda.em <-
function (documents, K, vocab, num.e.iterations, num.m.iterations,
    alpha, eta, annotations, params, variance, logistic = FALSE,
    lambda = 10, regularise = FALSE, method = "sLDA", trace = 0L, MaxNWts=3000)
{
    if (K > length(documents) && !regularise) {
        stop("K must be less than or equal to the length of documents.")
    }
    if (min(sapply(documents, length)) <= 0) {
        stop("All documents must have positive length.")
    }
    estimate.params <- function(document.sums, num.topics, logistic, MaxNWts=3000) {
        z.bar. <- t(document.sums)/colSums(document.sums)
        if (logistic) {
            if (!(is.integer(annotations)) || sum(0:(length(unique(annotations))-1) %in% unique(annotations)) != length(unique(annotations)))
            stop("Annotations must be consecutive integers starting from zero when sLDA and logistic.")
            model.fit <- nnet::multinom(annotations ~ z.bar. + 0, family = multinom(),MaxNWts=MaxNWts, trace=trace)
        }
        else if (regularise) {
            model.fit <- penalized::penalized(annotations ~ z.bar., unpenalized = ~0,
                lambda2 = 1/lambda^2, trace = FALSE)
        }
        else {
            model.fit <- lm(annotations ~ z.bar. + 0)
        }
        list(model = model.fit, coefs = as.double(t(coef(model.fit))))
    }

    method <- rep(method, length.out = num.m.iterations)
    variance <- rep(variance, length.out = num.m.iterations)
    num.e.iterations <- rep(num.e.iterations, length.out = num.m.iterations)
    result <- .slda.collapsed.gibbs.sampler(documents, K, vocab,
        num.e.iterations[1], alpha, eta, annotations, params,
        variance[1], logistic, method[1], lambda = lambda, trace = trace)
    fit <- estimate.params(result$document_sums, K, logistic,MaxNWts)
    params <- fit$coefs
    for (ii in seq(length.out = num.m.iterations - 1)) {
        result <- .slda.collapsed.gibbs.sampler(documents, K,
            vocab, num.e.iterations[ii + 1], alpha, eta, annotations,
            params, variance[ii + 1], logistic, method[ii + 1],
            lambda = lambda, initial = list(assignments = result$assignments),
            trace = trace)
        fit <- estimate.params(result$document_sums, K, logistic,MaxNWts)
        params <- fit$coefs
    }
    c(result, fit)
}
