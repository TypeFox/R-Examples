plotautocor <- function(data, ask=TRUE, lag.max=100, ...)
{
    if(!is.data.frame(data))
        data <- read.table(data, header=TRUE)

    data <- data[,-1]
    data.mcmc <- coda::as.mcmc(data)

    coda::autocorr.plot(data.mcmc, ask=ask, lag.max=lag.max, ...)     

    return(invisible())
}  
