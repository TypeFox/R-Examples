mixture <-
function (intensities, num.class, mix.method, threshold.0, threshold.k,
    mu.ini, sigma.ini, pi.ini, var.equal)
{
    method <- charmatch(mix.method, c("mixdist", "mclust", "EMmixt"))
    miss.threshold0 <- missing(threshold.0)
    miss.thresholdk <- missing(threshold.k)
    if (!miss.threshold0 && threshold.0 > min(intensities))
        num.class <- num.class - 1
    else {
        if (!miss.threshold0 && threshold.0 <= min(intensities))
            warning("threshold.0 ignored because it's smaller than intensity minimum value")
        threshold.0 <- -Inf
    }
    if (!miss.thresholdk && threshold.k < max(intensities))
        num.class <- num.class - 1
    else {
        if (!miss.thresholdk && threshold.k >= max(intensities))
            warning("threshold.k ignored because it's bigger than intensity maximum value")
        threshold.k <- Inf
    }
    yy <- intensities[intensities > threshold.0 & intensities < threshold.k]
    if (max(num.class) < 2) {
        cc <- unique(c(-Inf, threshold.0, threshold.k, Inf))
        out <- as.integer(cut(intensities, cc))
        attr(out, "mixture") <- NULL
        attr(out, "means") <- mean(yy)
        attr(out, "sds") <- 0
        attr(out, "probabilities") <- sapply(seq(along=unique(out)), function(j) ifelse(out==j,1,0))
        num.class <- 1
    }
    else {
        res <- mix(yy, method, num.class, mu.ini, sigma.ini, pi.ini, var.equal)
        num.class <- res$G
        out <- classCNV(intensities, res, threshold.0, threshold.k)
        attr(out, "mixture") <- res
        attr(out, "means") <- c(res$parameter$mu)
        attr(out, "sds") <- c(res$parameter$sigma)
    }
    attr(out, "meanRatio") <- intensities
    attr(out, "num.copies") <- c(1:num.class)
    if (!miss.threshold0) {
        attr(out, "means") <- c(min(intensities), attr(out, "means"))
        attr(out, "sds") <- c(0, attr(out, "sds"))
        attr(out, "num.copies") <- c(0, attr(out, "num.copies"))
    }
    if (!miss.thresholdk) {
        attr(out, "means") <- c(attr(out, "means"), max(intensities))
        attr(out, "sds") <- c(attr(out, "sds"), 0)
        attr(out, "num.copies") <- c(attr(out, "num.copies"), num.class + 1)
    }
    out
}

