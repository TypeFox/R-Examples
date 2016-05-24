brierscore <- function(object, data, group=NULL, decomp=FALSE, bounds=NULL, reverse=FALSE){

    ## Machinery to get group variable, if specified
    ## FIXME Use | within a formula for groups?
    if(missing(data)) data <- environment(object)
    args <- list(object=object, param=2, data=data, bounds=bounds,
                 reverse=reverse, decomp=decomp)

    if(!is.null(group)){
        mf <- match.call()
        m <- as.character(mf[match("group", names(mf), 0L)])
        group <- data[[m]]
        args <- c(args, list(group=group))
    }
    
    ## Get Brier scores + decomp using calcscore
    scores <- do.call("calcscore", args)

    ## If they supply group but not decomp, at least give them
    ## avg Briers
    if(decomp==FALSE & !is.null(group)){
        mnbrier <- tapply(scores, group, mean)
        scores <- list(rawscores=scores, mnbrier=mnbrier)
    }

    scores
}

bdecomp <- function(forecast, outcome, group){
    ## NB forecast is 2-column matrix.  First column is taken
    ## as outcome 1, second column is taken as outcome 0.
    f <- forecast[,1]
    d <- outcome==1

    mnbrier <- tapply((f-d)^2, group, mean)
    
    bias <- tapply(f, group, mean) - tapply(d, group, mean)

    slope <- as.numeric(diff(tapply(f, list(d, group), mean)))

    scatss <- tapply(f, list(d, group), function(x) (length(x)-1)*var(x))
    scatter <- 1/tapply(f, group, length) * apply(scatss,2, sum, na.rm=TRUE)

    data.frame(mnbrier, bias, slope, scatter)
}
