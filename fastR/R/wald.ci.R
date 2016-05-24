#' @export
wald.ci <-
function (x, n = 100, conf.level = 0.95) 
{
    alpha = 1 - conf.level
    p = x/n
    zstar <- -qnorm(alpha/2)
    interval <- p + c(-1, 1) * zstar * sqrt(p * (1 - p)/n)
    attr(interval, "conf.level") <- conf.level
    return(interval)
}
