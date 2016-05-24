predict.bma <-
function (object, newdata = NULL, exact = FALSE, topmodels = NULL, 
    ...) 
{
    if (!is.bma(object)) {
        stop("you need to provide a BMA object")
        return()
    }
    if (!is.null(topmodels)) {
        if (!(is.numeric(topmodels) && is.vector(topmodels))) {
            stop("topmodels must denote the models to take into account, e.g. 1:5 for the best five.")
        }
        else if (object$topmod$nbmodels < max(topmodels)) {
            stop(paste("Only", object$topmod$nbmodels, "best models are available, but you asked to take the", 
                max(topmodels), "-best model into account."))
        }
        object = object[unique(topmodels)]
    }
    if ((!missing(topmodels)) && missing(exact)) 
        exact = TRUE
    betas = estimates.bma(object, exact = exact, order.by.pip = FALSE, 
        include.constant = FALSE, std.coefs = FALSE, condi.coef = FALSE)[, 
        2]
    if (is.null(newdata)) {
        newX <- as.matrix(object$X.data[, -1, drop = FALSE])
    }
    else {
        newX = as.matrix(newdata)
        if (!is.numeric(newX)) 
            stop("newdata must be numeric!")
        if (is.vector(newdata)) 
            newX = matrix(newdata, 1)
        if (ncol(newX) != length(betas)) {
            if (ncol(newX) == length(betas) + 1) {
                newX = newX[, -1, drop = FALSE]
            }
            else {
                stop("newdata must be a matrix or data.frame with", 
                  length(betas), "columns.")
            }
        }
        orinames = colnames(object$X.data[, -1, drop = FALSE])
        if (!is.null(colnames(newX)) && !is.null(orinames)) {
            if (all(orinames %in% colnames(newX)) && !all(orinames == 
                colnames(newX))) {
                warning("argument newdata had to be reordered according to its column names. Consider submitting the columns of newdata in the right order.")
                newX = newX[, orinames, drop = FALSE]
            }
        }
    }
    cons = .post.constant(object$X.data, betas)
    return(as.vector(newX %*% betas) + cons)
}
