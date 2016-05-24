nonzero.customizedGlmnet <-
function(object, lambda = NULL, ...)
{
    groups = as.character(sort(unique(object$groupid)))
    selected = list()

    for (group in groups) {

        selected[[group]] = predict(object$fit[[group]],
            s = lambda/object$fit[[group]]$nobs, type = "nonzero")
    }

    selected
}
