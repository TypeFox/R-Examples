#
#  Copyright (C) 2005-2012 Friedrich Leisch
#  $Id: stepflexclust.R 13 2013-07-01 05:15:35Z leisch $
#


stepFlexclust <- function(x, k, nrep=3, verbose=TRUE,
                          FUN=kcca, drop=TRUE, group=NULL,
                          simple=FALSE, save.data=FALSE, seed=NULL,
                          multicore=TRUE, ...)
{
    MYCALL <- match.call()

    if(is.character(FUN)) FUN <- get(FUN, mode="function")
    
    if(!is.null(seed)) set.seed(seed)
    
    if(!is.logical(multicore))
        stop("argument multicore is not logical (TRUE or FALSE)")
    
    bestKcca <- function(x, k, ...)
    {
        seed <- as.list(round(2^31 * runif(nrep, -1, 1)))

        res <- MClapply(seed, 
                        function(y){
                            set.seed(y)
                            if(verbose) cat(" *")
                            FUN(x=x, k=k, group=group, simple=TRUE,
                                save.data=FALSE,
                                ...)
                        }, multicore=multicore)

        distsum <- sapply(res, function(y) info(y, "distsum"))
        res[[which.min(distsum)]]
    }
    
    k = as.integer(k)
    if(length(k)==0)
        return(list())
    
    z = list()
    MYCALL1 <- MYCALL
    if("drop" %in% names(MYCALL))
        MYCALL1[["drop"]] <- NULL
    
    for(n in 1:length(k)){
        if(verbose) cat(k[n], ":")
        kn <- as.character(k[n])
        z[[kn]] = bestKcca(x=x, k=k[n], ...)
        MYCALL1[["k"]] <- k[n]
        z[[kn]]@call <- MYCALL1

        if(!simple){
            ## x is usually at the beginning of kcca() pre-porcessed,
            ## here we have to do it manually!
            z[[kn]] <- simple2kcca(x=z[[kn]]@family@preproc(x),
                                   from=z[[kn]], group=group)
        }
        else{
        }
        if(verbose) cat("\n")
    }

    if(save.data){
        me <- ModelEnvMatrix(designMatrix=x)
        for(n in seq_along(z))
            z[[n]]@data <- me
    }
    
    if(drop && length(k)==1){
        return(z[[1]])
    }
    else{
        z <- new("stepFlexclust", models=z, k=as(k, "integer"),
                 nrep=as(nrep, "integer"), call=MYCALL)
        if(simple){
            x <- z@models[[1]]@family@preproc(x)
            z@xcent <- z@models[[1]]@family@cent(x)
            z@totaldist <-
                sum(z@models[[1]]@family@dist(x,
                                              matrix(z@xcent,nrow=1)))
        }
        else{
            z@xcent <- z@models[[1]]@xcent
            z@totaldist <- z@models[[1]]@totaldist
        }
        return(z)
    }
}

stepcclust <- function(...) stepFlexclust(..., FUN=cclust)


setMethod("show", "stepFlexclust",
function(object)
{
    cat("stepFlexclust object of family",
        sQuote(object@models[[1]]@family@name),"\n\n")
    cat("call:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\n")
    
    z <- data.frame(iter = sapply(object@models, function(x) x@iter),
                    converged = sapply(object@models, function(x) x@converged),
                    distsum = sapply(object@models,
                                     function(x) info(x, "distsum")))

    z1 <- data.frame(iter = NA,
                     converged = NA,
                     distsum = object@totaldist)
    
    z <- rbind(z1, z)
    
    print(z, na.string="")
})
    
setMethod("getModel", "stepFlexclust",
function(object, which=1)
{
    object@models[[which]]
})

setMethod("[[", signature(x="stepFlexclust", i="ANY", j="missing"),
function(x, i, j) getModel(x, i))

###**********************************************************

setGeneric("clusterJaccard",
function(object, object2, ...) standardGeneric("clusterJaccard"))

setMethod("clusterJaccard",
          signature(object="kccasimple", object2="kccasimple"),
function(object, object2, ...)
clusterJaccard(object@cluster, object2@cluster))

setMethod("clusterJaccard",
          signature(object="integer", object2="integer"),
function(object, object2, ...)
{
    k1 <- max(object)
    k2 <- max(object2)
    z <- matrix(double(1), nrow=k1, ncol=k2)
    for(m in 1:k1){
        ok1 <- (object==m)
        for(n in 1:k2){
            ok2 <- (object2==n)
            z[m,n] <- sum(ok1&ok2)/sum(ok1|ok2)
        }
    }
    z
})

