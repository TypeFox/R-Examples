convert.mcmc.list <- function(x){
    if (!is.mcmc.list(x)){
        if (!is.mcmc(x)){
            if ("list" %in% class(x)){
                x <- lapply(x, as.mcmc)
            } else {
                x <- as.mcmc(x)
            }
        }
        if (!is.mcmc.list(x)){
            x <- mcmc.list(x)
        }
    }
    return(x)
}
