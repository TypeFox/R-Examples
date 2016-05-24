predict.ewkm <- function (object, data, ...)
{
    cluster.names <- rownames(object$centers)
    cluster.vars <- colnames(object$centers)
    out <- apply(data[cluster.vars], 1,
                 function(d) cluster.names[which.min(lapply(1:nrow(object$centers),
                                                            function(i) sqrt(sum(object$weights[i,] * abs(d - object$centers[i,])^2))))])
    out <- sapply(out, function(x) ifelse(length(x), x, NA))
    return(out)
}
