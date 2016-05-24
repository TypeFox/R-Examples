sphscore <- function(object, data, group=NULL, bounds=NULL, reverse=FALSE){

    ## Machinery to get group variable, if specified
    ## FIXME Use | within a formula for groups?
    if(missing(data)) data <- environment(object)
    args <- list(object=object, param=2, fam="sph", data=data,
                 bounds=bounds, reverse=reverse)

    if(!is.null(group)){
        mf <- match.call()
        m <- as.character(mf[match("group", names(mf), 0L)])
        group <- data[[m]]
        args <- c(args, list(group=group))
    }
    
    ## Get sph scores + decomp using calcscore
    scores <- do.call("calcscore", args)

    if(!is.null(group)){
        mnsph <- tapply(scores, group, mean)
        scores <- list(rawscores=scores, mnsph=mnsph)
    }

    scores
}
