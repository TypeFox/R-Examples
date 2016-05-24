`predict.prcurve` <- function(object, newdata, ...) {
    if(missing(newdata))
        return(fitted(object))

    ## check the variable names in newdata match with original data
    ## essentially, this will only work if the names match, hence
    ## join() likely useful for the user
    nNew <- colnames(newdata)
    nData <- colnames(object$data)
    if (!isTRUE(all.equal(nNew, nData))) {
        if (isTRUE(all.equal(sort(nNew), sort(nData)))) {
            newdata <- newdata[, nData]
        } else {
            stop("Variables in 'newdata' don't match with training datat.")
        }
    }

    ## otherwise project points on to the curve
    p <- get.lam(data.matrix(newdata), s = object$s, tag = object$tag,
                 stretch = object$stretch)
    out <- p$s
    attr(out, "tag") <- p$tag
    out
}
