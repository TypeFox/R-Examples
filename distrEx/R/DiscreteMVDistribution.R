# are rows of x also rows of y
# if yes, row number is computed
.rowNr <- function(x, y){
    x <- as.matrix(x)
    y <- as.matrix(y)
    if(ncol(x) != ncol(y))
        stop("'x' and 'y' have different numbers of columns")
    res <- numeric(0)
    DistrResolution <- getdistrOption("DistrResolution")
    for(i in 1:nrow(x))
        for(j in 1:nrow(y))
            if(isTRUE(all.equal(x[i,], y[j,], tolerance = DistrResolution))) 
               res <- c(res, j)
    return(res)
}

# Generating function
DiscreteMVDistribution <- function(supp, prob, Symmetry = NoSymmetry()){
    if(!is.numeric(supp)) 
        stop("'supp' has to be numeric")
    if(!is.matrix(supp)){
        supp <- t(as.matrix(supp))
    }
    if(any(!is.finite(supp)))
        stop("inifinite or missing values in supp")
    len <- nrow(supp)
    if(missing(prob)){
        prob <- rep(1/len, len)
    }else{
        if(len != length(prob))
            stop("number of columns of 'supp' != length of 'prob'")
        if(any(!is.finite(prob)))
            stop("inifinite or missing values in prob")
        if(sum(prob) != 1)
            stop("sum of 'prob' != 1")
        if(!all(prob >= 0))
            stop("'prob' contains values < 0")
    }
    if(any(duplicated(supp))){
        warning("collapsing to unique support values")
        usupp <- supp[!duplicated(supp),]
        ind <- .rowNr(supp, usupp)
        prob <- as.vector(tapply(prob, ind, sum))
        supp <- usupp
        len <- nrow(supp)
    }

    rfun <- function(n){ 
        ind <- sample(x = 1:len, size = n, replace = TRUE, prob = prob)
        return(supp[ind,])
    }
  
    dfun <- function(x){ 
        if(is.vector(x)) x <- t(x)
        ind <- .rowNr(x, supp)
        res <- numeric(nrow(x))
        if(length(ind) == 0)
            return(res)
        else{
            p <- prob[ind]
            if(length(p) < nrow(x)){
                res[ind] <- prob[ind]
                return(res)
            }else return(p)
        }
    }
  
    pfun <- function(lower, upper){ 
        if(!is.numeric(lower) || !is.numeric(upper))
            stop("'lower' and 'upper' have to be numeric vectors")
        if(length(lower) != ncol(supp))
            stop("wrong dimension of 'lower'")
        if(length(upper) != ncol(supp))
            stop("wrong dimension of 'upper'")
        ind1 <- apply(t(supp) >= lower, 2, all)
        ind2 <- apply(t(supp) <= upper, 2, all)
        ind <- ind1 & ind2
        sum(prob[ind])
    }
        
    MVD <- new("DiscreteMVDistribution")
    MVD@r <- rfun
    MVD@d <- dfun
    MVD@p <- pfun
    MVD@q <- NULL
    MVD@param <- NULL
    MVD@img <- EuclideanSpace(dimension = floor(ncol(supp)))
    MVD@support <- supp
    MVD@.withSim <- FALSE 
    MVD@.withArith <- FALSE
    MVD@.logExact <- TRUE 
    MVD@.lowerExact <- FALSE
    MVD@Symmetry <- Symmetry
    
    return(MVD)
}

setMethod("support", "DiscreteMVDistribution", function(object) object@support)
setMethod("dim", "DiscreteMVDistribution", function(x)ncol(x@support))
