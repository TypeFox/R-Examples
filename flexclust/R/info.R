#
#  Copyright (C) 2005 Friedrich Leisch
#  $Id: info.R 3 2013-06-12 10:06:43Z leisch $
#

setMethod("info", signature(object="flexclust", which="character"),
function(object, which, drop=TRUE, ...)
{
    INFOS <- c(names(object@clusinfo))
    if(nrow(object@cldist) | all(c("size","av_dist") %in% INFOS))
        INFOS <- c(INFOS, "distsum")
    
    if("help" %in% which){
        return(INFOS)
    }

    which <- INFOS[pmatch(which, INFOS)]
    
    if(any(is.na(which))){
        stop(paste("Requested info not available, use",
                   sQuote(paste("which=",dQuote("help"),sep="")),
                   "to list available infos."))
    }

    if(any(which=="distsum")){
        if(all(c("size","av_dist") %in% INFOS))
            return(sum(object@clusinfo$size *
                       object@clusinfo$av_dist))
        else
            return(sum(object@cldist[,1]))
    }
    else{
        z <- object@clusinfo[,which,drop=drop]
        if(is.vector(z))
            names(z) <- rownames(object@clusinfo)
        return(z)
    }
})

###**********************************************************

setMethod("parameters", signature(object="kccasimple"),
function(object, ...)
{
    object@centers
})


###**********************************************************

setGeneric("clusterSim",
           function(object, ...) standardGeneric("clusterSim"))

setMethod("clusterSim", signature(object="kcca"),
function(object, data=NULL, method=c("shadow", "centers"),
         symmetric=FALSE, ...)
{
    method <- match.arg(method)

    if((method=="shadow") && is.null(data)){
        z <- object@clsim
        if(symmetric) z <- (z+t(z))/2
    }
    else{
        z <- callNextMethod(object=object, data=data, method=method,
                            symmetric=symmetric, ...)
    }

    z
})

setMethod("clusterSim", signature(object="kccasimple"),
function(object, data=NULL, method=c("shadow", "centers"),
         symmetric=FALSE, ...)
{
    method <- match.arg(method)
    
    if(object@k==1) return(matrix(1))

    if((method=="shadow")){
        if(is.null(data)) data <- getData(object)

        if(any(is.na(object@cluster)))
            data <- data[!is.na(object@cluster),]

        distmat <- object@family@dist(data, object@centers)
        cluster <- object@family@cluster(n=2, distmat=distmat)
        z <- flexclust:::computeClusterSim(distmat, cluster)
        if(symmetric) z <- (z+t(z))/2
    }
    else{
        z <- object@family@dist(object@centers, object@centers)
        z <- 1-z/max(z)
    }
        
    z
})
