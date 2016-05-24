gx.scores <-
function (xx, tholds, rwts = NULL, setna = FALSE) 
{
    if (!is.matrix(xx)) 
        stop(deparse(substitute(xx)), " is not a Matrix")
    temp.x <- remove.na(xx)
    x <- temp.x$x
    n <- temp.x$n
    ncolx <- temp.x$m
    nthld <- length(tholds)
    if (ncolx != nthld) 
        stop("\n  Number of variables and thresholds do not match")
    if (!is.null(rwts) && length(rwts) != nthld)
        stop("\n  Number of thresholds and weights do not match")
    for (i in 1:n) {
        for (j in 1:nthld) {
            x[i,j] <- x[i,j]/tholds[j]
            if (x[i,j] < 1) x[i,j] <- 0
            if (!is.null(rwts)) x[i,j] <- x[i,j] * rwts[j]
        }
    }
    scores <- rowSums(x)
    scores[scores <= 0] <- 0
    if (setna) scores[scores <= 0] <- NA
    invisible(list(input = deparse(substitute(xx)), tholds = tholds, 
        rwts = rwts, scores = scores))
}
