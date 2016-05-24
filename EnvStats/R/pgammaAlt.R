pgammaAlt <-
function (q, mean, cv = 1, lower.tail = TRUE, log.p = FALSE) 
{
    shape <- cv^-2
    scale <- mean/shape
    stats::pgamma(q = q, shape = shape, scale = scale, lower.tail = lower.tail, 
        log.p = log.p)
}
