#
#  Copyright (C) 2005-2008 Friedrich Leisch
#  $Id: utils.R 3 2013-06-12 10:06:43Z leisch $
#

list2object = function(from, to){
    n = names(from)
    s = slotNames(to)
    p = pmatch(n, s)
    if(any(is.na(p)))
        stop(paste("\nInvalid slot name(s) for class",
                   to, ":", paste(n[is.na(p)], collapse=" ")))
    names(from) = s[p]
    do.call("new", c(from, Class=to))
}

printIter <- function(iter, logLik, label="Log-likelihood",
                      format="f", width=12)
    cat(formatC(iter, width=6),
        label, ":", formatC(logLik, width=width, format=format),"\n")
    



## library(colorspace)
## ORDER=c(1,3,5,7,2,4,6,8)
## dput(x[ORDER])

## x = hcl(seq(0, 360*7/8, length = 8), c=30, l=85)
LightColors <- c("#FAC8D1", "#D4D8AE", "#A3E0D8", "#D5D0F6",
                 "#EECEB7", "#B5DFBD", "#B2DAEF", "#F1C8EA")

## dput(hcl(seq(0, 360*7/8, length = 8), c=65, l=85)[ORDER])
MedColors <- c("#FFB8CC", "#D4DB76", "#2BEDDC", "#D5CBFF",
               "#FFC88F", "#88E99F", "#72E2FF", "#FFB7FF")


## x = hcl(seq(0, 360*7/8, length = 8), c=100, l=65)
FullColors <- c("#FF6C91", "#9DA700", "#00C1A9", "#9F8CFF",
                "#DE8C00", "#00BA38", "#00B4F0", "#F564E3")

##  x=hcl(seq(0, 360*7/8, length = 8), c=40, l=65)
DarkColors <- c("#CC8D99", "#9DA268", "#4EADA2", "#9E98CA",
                "#BE9675", "#71AB7E", "#69A6C0", "#C28DBA")




flxColors <- function(n=1:8, color=c("full","medium", "light","dark"),
                      grey=FALSE)
{
    color <- match.arg(color)
    
    if(color=="light"){
        if(grey)
            return("#D4D4D4")
        else
            return(LightColors[n])
    }
    if(color=="medium"){
        if(grey)
            return("#D4D4D4")
        else
            return(MedColors[n])
    }
    else{
        if(grey) return("#9E9E9E")
        
        if(color=="full"){
            return(FullColors[n])
        }
        else{
            return(DarkColors[n])
        }
    }
}

###**********************************************************

getData <- function(x, error=FALSE)
{
    if(empty(x@data)){
        if(error) stop("Cluster object contains no data.")
        z <- NULL
    }
    else{
        z <- x@data@get("designMatrix")
    }
    z
}

###**********************************************************

## if length(col)<=k first recycle to k, then do col[cluster]
## else simply recycle to number of observations
expandColors <- function(col, object)
{
    k <- object@k
    
    if(is.null(col))
        col <- flxColors(n=1:min(k, 8) , color="full")
    
    if(length(col) <= k){
        col <- rep(col, length=k)
        col <- col[object@cluster]
    }
    else{
        col <- rep(col, length=nrow(object@cldist))
    }
    
    col
}

###**********************************************************

MClapply <- function(X, FUN, multicore=TRUE, ...)
{
    if(inherits(multicore, "cluster"))
        parLapply(multicore, X, FUN)
    else if(multicore)
        mclapply(X, FUN, ...)
    else
        lapply(X, FUN, ...)
}

