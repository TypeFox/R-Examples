predict.customizedGlmnet <-
function(object, lambda, ...)
{
    groups = as.character(sort(unique(object$groupid)))
    prediction = matrix(NA, nrow(object$x$test), length(lambda))

    for (group in groups) {

        x = object$x$test[object$groupid == group, ]

        if (sum(object$groupid == group) == 1) {
            x = t(x)
        }

        type = "response"
        if (is.element(object$family, c("binomial", "multinomial"))) {
            type = "class"
        }

        prediction[object$groupid == group, ] = predict(object$fit[[group]],
            x, s = lambda/object$fit[[group]]$nobs, type = type)
    }

    prediction
}
