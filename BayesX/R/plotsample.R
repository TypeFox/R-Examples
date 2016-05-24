plotsample <- function(data, ask=TRUE, ...)
{
    if(!is.data.frame(data))
        data <- read.table(data, header=TRUE)
    
    old.par <- par(no.readonly=TRUE)	
    on.exit(par(old.par))
    
    i <- 1
    count <- 0
    while(i < length(data)){
        if(count %% 6 == 0){
            count <- 0
            par(mfrow=c(3, 2), omi=c(1.5, 1, 0, 1.5)/2.54)
        }
        plot(data[, 1], data[, i + 1], type="l", xlab="iteration", ylab=dimnames(data)[[2]][i + 1], ...)
        par(ask=ask)
        count <- count + 1
        i <- i + 1
    }


    return(invisible())
}  

plotsample.coda <- function(data, ask=TRUE, ...)
{
    if(!is.data.frame(data))
        data <- read.table(data, header=TRUE)

    data <- data[,-1]
    data.mcmc <- coda::as.mcmc(data)
    plot(data.mcmc, ask=ask, ...)     

    return(invisible())
}  
