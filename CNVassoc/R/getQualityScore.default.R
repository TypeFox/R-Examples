getQualityScore.default <-
function (x, sds, w, type, iter = 10000, threshold = 0.1, ...)
{
    mu <- x
    J <- length(mu)
    qs.type <- charmatch(tolower(type), tolower(c("class", "CNVtools", "CANARY")))
    if (is.na(qs.type))
        stop(" argument 'type' must be either 'class' , 'CNVtools' or 'CANARY'")
    if (qs.type == 1) {
        p <- c()
        for (j in 1:J) {
            X <- rnorm(iter, mu[j], sds[j])
            Y <- sapply(1:J, function(s) w[s] * dnorm(X, mu[s], sds[s]))
            p <- c(p, mean(apply(Y, 1, which.max) == j))
        }
        out <- sum(p * w)
    }
    if (qs.type == 2) {
        sds <- sds[order(mu)]
        w <- w[order(mu)]
        mu <- sort(mu)
        dmu <- abs(mu[1:(J - 1)] - mu[2:J])
        av.sds <- (w[1:(J - 1)] * sds[1:(J - 1)] + w[2:J] * sds[2:J])/(w[1:(J - 1)] + w[2:J])
        weights <- w[1:(J - 1)] * w[2:J]
        out <- sum(weights * dmu/av.sds)/sum(weights)
    }
    if (qs.type == 3) {
        f <- function(x) {
            Y <- sapply(1:J, function(s) w[s] * dnorm(x, mu[s], sds[s]))
            index <- which.max(Y)
            max1 <- Y[index]
            max2 <- max(Y[-index])
            ratio <- max2/max1
            sum(Y) * as.integer(ratio > threshold)
        }
        minim <- which.min(mu)
        maxim <- which.max(mu)
        limits <- c(mu[minim] - 3 * sds[minim], mu[maxim] + 3 * sds[maxim])
        fapply <- function(x) sapply(x, f)
        out <- integrate(fapply, limits[1], limits[2])$value
        attr(out, "threshold") <- threshold
    }
    attr(out, "type") <- qs.type
    return(out)
}
