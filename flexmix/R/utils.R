#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: utils.R 4878 2013-02-08 09:09:35Z gruen $
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

printIter = function(iter, logLik, label="Log-likelihood")
    cat(formatC(iter, width=4),
        label, ":", formatC(logLik, width=12, format="f"),"\n")
    

## library(colorspace)
## dput(x[c(1,3,5,7,2,4,6,8)])

## x = hcl(seq(0, 360*7/8, length.out = 8), c=30)
LightColors <- c("#F9C3CD", "#D0D4A8", "#9DDDD5", "#D1CCF5",
                 "#EDCAB2", "#AFDCB8", "#ACD7ED", "#EFC4E8")
    
## x = hcl(seq(0, 360*7/8, length.out = 8), c=100, l=65)
FullColors <- c("#FF648A", "#96A100", "#00BCA3", "#9885FF",
                "#DC8400", "#00B430", "#00AEEF", "#F45BE1")


###**********************************************************

## similar defaults to silhouette plots in flexclust
unipolarCols <- function(n, hue=0, chr=50, lum = c(55, 90))
{
    lum <- seq(lum[1], lum[2], length=n)
    hcl(hue, chr, lum)
}


bipolarCols <- function(n, hue=c(10, 130), ...)
{        
    if(n%%2){ # n odd
        n2 <- (n-1)/2
        c1 <- unipolarCols(n2, hue[1])
        c2 <- rev(unipolarCols(n2, hue[2]))
        return(c(c1, "white", c2))
    }
    else{ # n even
        n2 <- n/2
        c1 <- unipolarCols(n2, hue[1])
        c2 <- rev(unipolarCols(n2, hue[2]))
        return(c(c1, c2))
    }
}

###**********************************************************
