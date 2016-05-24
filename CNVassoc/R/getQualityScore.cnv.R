getQualityScore.cnv <-
function (x, type = "class", iter = 10000, threshold = 0.1, ...)
{
    obj <- x
    if (is.null(attr(obj, "means"))) {
        ff <- function(threshold) {
            res <- apply(attr(obj, "probabilities"), 1, function(x) {
                index <- which.max(x)
                max1 <- x[index]
                max2 <- max(x[-index])
                (max2/max1) > threshold
            })
            mean(res)
        }
        out <- ff(threshold)
        attr(out, "type") <- 0
        attr(out, "threshold") <- threshold
    }
    else {
        mu <- attr(obj, "means")
        sds <- attr(obj, "sds")
        w <- attr(obj, "pi")
        batches <- attr(x, "batches")
        if (is.null(batches))
          out <- getQualityScore(mu, sds, w, type, iter = iter, threshold = threshold)
        else{
          out <- list()
          bb <- sort(unique(batches))
          for (i in bb)
            out[[i]] <- getQualityScore(mu[i,], sds[i,], w[i,], type, iter = iter, threshold = threshold)
          attr(out,"batch") <- bb
          attr(out, "type") <- attr(out[[1]], "type")
          attr(out, "threshold") <- threshold
        }
    }
    class(out) <- "QualityScore"
    return(out)
}
