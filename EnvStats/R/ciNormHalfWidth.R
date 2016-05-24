ciNormHalfWidth <-
function (n.or.n1, n2 = n.or.n1, sigma.hat = 1, conf.level = 0.95, 
    sample.type = ifelse(missing(n2), "one.sample", "two.sample")) 
{
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    if (!is.vector(n.or.n1, mode = "numeric") || !is.vector(sigma.hat, 
        mode = "numeric") || !is.vector(conf.level, mode = "numeric")) 
        stop("'n.or.n1', 'sigma.hat', and 'conf.level' must be numeric vectors.")
    if (!all(is.finite(n.or.n1)) || !all(is.finite(sigma.hat)) || 
        !all(is.finite(conf.level))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'n.or.n1', 'sigma.hat', or 'conf.level'"))
    if (any(n.or.n1 < 2)) 
        stop("All values of 'n.or.n1' must be greater than or equal to 2")
    if (any(sigma.hat < .Machine$double.eps)) 
        stop("All values of 'sigma.hat' must be positive.")
    if (any(conf.level <= .Machine$double.eps) || any(conf.level >= 
        1 - .Machine$double.eps)) 
        stop("All values of 'conf.level' must be greater than 0 and less than 1")
    if (sample.type == "two.sample") {
        if (!missing(n2)) {
            if (!is.vector(n2, mode = "numeric")) 
                stop("'n2' must be a numeric vector")
            if (!all(is.finite(n2))) 
                stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
                  "Undefined (Nan) values are not allowed in 'n2'"))
            if (any(n2 < 2)) 
                stop("All values of 'n2' must be greater than or equal to 2")
        }
        half.width <- qt(1 - (1 - conf.level)/2, n.or.n1 + n2 - 
            2) * sigma.hat * sqrt(1/n.or.n1 + 1/n2)
    }
    else half.width <- (qt(1 - (1 - conf.level)/2, n.or.n1 - 
        1) * sigma.hat)/sqrt(n.or.n1)
    names(half.width) <- NULL
    half.width
}
