hlnorm <- function(x, meanlog = 0, sdlog = 1,
                   shape = 1 / sdlog, scale = exp(meanlog), prop = 1,
                   log = FALSE){
    ## 'prop' added in version 2.3-0. 
    ## shape = 1 / sdlog, scale = exp(meanlog)
    meanlog <- log(scale)
    sdlog <- 1 / shape
    ret <- prop * dlnorm(x, meanlog, sdlog) /
        plnorm(x, meanlog, sdlog, lower.tail = FALSE)
    if (log) ret <- ifelse(ret <= 0, -Inf, log(ret))
    return (ret)

}

Hlnorm <- function(x, meanlog = 0, sdlog = 1,
                   shape = 1 / sdlog, scale = exp(meanlog), prop = 1,
                   log.p = FALSE){
    ## 'prop' added in version 2.3-0. 
    ## shape = 1 / sdlog, scale = exp(meanlog)
    meanlog <- log(scale)
    sdlog <- 1 / shape
    ret <- -prop * plnorm(x, meanlog, sdlog, lower.tail = FALSE, log.p = TRUE)
    if (log.p) ret <- log(ret)
    return (ret)
}
