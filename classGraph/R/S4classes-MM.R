### This was
## ("/u/maechler/R/Meetings-Kurse-etc/2005-DSC/ S4classes.R")
## BioC (Robert Gentleman): class2Graph() and utilities
##
## but that -- and also the 'graph' package -- had bugs!


###--- 2nd, the 'S4classes' utilites corrections : -------------------------

fullyQcName <- function(x) {
    pName <- attr(x, "package")
    if (is.null(pName)) x else paste(pName, x, sep = ":")
}

superClasses <- function(x) {
    if(!is(x, "classRepresentation") )
        return("must have a class representation object")
    superCs  <- names(x@contains)
    if(length(superCs) == 0 )
        return(character(0))
    directSCs  <- sapply(x@contains,
                         function(x) if(length(x@by) > 0 ) FALSE else TRUE)
    pkgNames  <- sapply(x@contains, function(x) x@package)
    clss  <- superCs[directSCs]
    pkgNames  <- pkgNames[directSCs]
    ans  <- vector("list", length = length(clss))
    for( i in 1:length(clss)) {
        v  <- clss[i]
        attr(v, "package") <- pkgNames[i]
        ans[[i]] <- v
    }
    return(ans)
}

### FIXME: this must have a bug too,
### ----- since (cg2 <- class2Graph("dtrMatrix", fullNames = FALSE))
### is almost empty;

### No, actually, the culprit is
## >>  getAllSuperClasses(getClass("dtrMatrix"))
## which returns an empty character vector
## even though dtrMatrix does have several superclasses;
## namely "dgeMatrix" `` directly, with explicit coerce ''
## {and 4 more via "dgeMatrix"} : but actually

## MM: use 'package' and  getClass(*, where=.) ! to find private classes
class2Graph <-
    function(class, fullNames = TRUE, simpleOnly = FALSE, bottomUp = FALSE,
             package = class@package)
{
    if(is(class, "character"))
	class <- getClass(class)
    if( !is(class, "classRepresentation") )
        stop("need a character or a classRepresentation")

    cname  <- as.character(class@className)
    where <- asNamespace(package)
    superClasses  <- getAllSuperClasses(class, simpleOnly = simpleOnly)
    ## MM:                                      ^^^^^^^^^^^^^^^^^^^^^ important

    ##handle the one node graph separately
    if( length(superClasses) == 0 ) {
        eL  <- setNames(list(numeric(0)), cname)
        return(new("graphNEL", edgeL = eL, nodes = cname))
    }
    ##otherwise build a simple incidence matrix
    nN  <- length(superClasses)+1
    rmat  <- matrix(0, nN, nN)
    dimnames(rmat) <-
        list(c(cname, superClasses),
             c(cname, superClasses))
    sCn  <- superClasses(class)
    fNms  <- rep("", nN)
    if( fullNames )
        fNms[1]  <- fullyQcName(class@className)
    rmat[cname, as.character(sCn)] <- 1
    for(i in 1:(nN-1)) {
        tc  <- getClass(superClasses[i], where=where)
        tCn  <- superClasses(tc)
        rmat[superClasses[i], as.character(tCn)] <- 1
        if(fullNames)
            fNms[i+1] <- fullyQcName(tc@className)
    }
    if (fullNames)
        dimnames(rmat) <- list(fNms, fNms)
    return(as(if(bottomUp) t(rmat) else rmat, "graphNEL"))
}
