`residuals.prcurve` <- function(object,
                                which = c("distance","raw","smooths","pca"),
                                ...) {
    which <- match.arg(which)

    ## predict locations of curve in PCA space and exctract the
    ## site scores for each sample
    if( isTRUE(all.equal(which, "pca"))) {
        p <- predict(object[["ordination"]], object[["s"]],
                     which = "wa", scaling = 1) ## site scaling
        scrs <- scores(object[["ordination"]], scaling = 1,
                       choices = seq_along(eigenvals(object[["ordination"]])),
                       display = "sites")
    }

    res <- switch(which,
                  distance = diag(distance(object$s, object$data,
                  method = "euclidean")),
                  raw = object[["data"]] - object[["s"]],
                  smooths = sapply(object$smooths, function(x, ...)
                  residuals(x$model, ...), ...),
                  pca = p - scrs)
    res
}
