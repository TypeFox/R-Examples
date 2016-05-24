#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* campus *.* lmu *.* de]
## Project: BayesX
## Time-stamp: <[nbAndGraConversion.R] by DSB Son 22/02/2009 19:10 (CET) on daniel@puc-home>
##
## Description:
## Convert an nb object (neighborhood structure) from package spdep to the graph format
## required by BayesX, and vice versa.
##
## History:
## 19/02/2009   file creation: first version (without optional weights, as weights
##              are not considered in read/write.gra)
## 22/02/2009   remove spdep dependence 
#####################################################################################


nb2gra <- function(nbObject)
{
    ## check if S3 class of nbObject is "nb"
    stopifnot(inherits(x=nbObject,
                       what="nb"))

    ## convert to (negative) binary neighbors matrix
    regionIds <- attr(nbObject, "region.id")
    ret <- matrix(data=0,
                  nrow=length(regionIds),
                  ncol=length(regionIds),
                  dimnames=
                  list(regionIds,
                       regionIds))

    for(i in seq_along(nbObject))
    {
        ret[i, nbObject[[i]]] <- - 1
    }


    ## and to gra format
    diag(ret) <- - rowSums(ret)
    class(ret) <- "gra"

    return(ret)
}


gra2nb <- function(graObject)
{
    ## check if S3 class of nbObject is "gra"
    stopifnot(inherits(x=graObject,
                       what="gra"))

    ## save region names and delete them
    ## (so that the list below will not have names attached)
    regionNames <- rownames(graObject)
    dimnames(graObject) <- NULL
    
    ## make list of neighbors
    ret <- apply(graObject,
                 MARGIN=1,
                 FUN=function(row) which(row == -1))

    ## attach necessary attributes
    ret <- structure(ret,
                     class="nb",
                     region.id=regionNames,
                     call=match.call(),
                     type="queen",
                     sym=TRUE)

    ## and return the nb object
    return(ret)
}



