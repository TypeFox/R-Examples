orlm.forboot.fixed <- function (data, indices, ...) 
{
    e <- data$e[indices]
    dat <- cbind(data$fit+e, data[,2:(ncol(data)-2)])
    aus <- orlm(lm(dat), ...)$b.restr
    return(aus)
}
