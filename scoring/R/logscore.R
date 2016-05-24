logscore <- function(object, data, group=NULL){

    ## Machinery to get group variable, if specified
    ## FIXME Use | within a formula for groups?
    if(missing(data)) data <- environment(object)
    args <- list(object=object, fam="log", data=data)

    if(!is.null(group)){
        mf <- match.call()
        m <- as.character(mf[match("group", names(mf), 0L)])
        group <- data[[m]]
        args <- c(args, list(group=group))
    }
    
    ## Get log scores + decomp using calcscore
    scores <- do.call("calcscore", args)

    if(!is.null(group)){
        mnlog <- tapply(scores, group, mean)
        scores <- list(rawscores=scores, mnlog=mnlog)
    }

    scores
}
