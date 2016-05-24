trimImage <- function(x, max2trim=.Machine$double.eps, na.rm=TRUE,
          returnIndices2Keep=FALSE, ...){
##
## 1.  check arguments
##
    if(!is.numeric(x)){
        stop('x must be numeric;  is ', class(x))
    }
#  1.1.  dim(x)?
    if(length(dim(x))<2){
        stop('length(dim(x)) = ', length(dim(x)),
             ' < 2; must be a matrix or 3-d array')
    }
    if(length(dim(x))>3){
        stop('length(dim(x)) = ', length(dim(x)),
             ' > 3;  must be a matrix or a 3-d array')
    }
#  1.2.  is.logical(na.rm)
    if(!is.logical(na.rm)){
        stop('na.rm must be logical;  class(na.rm) = ',
             class(na.rm))
    }
##
## 2.  is.list(returnIndices2Keep)?
##
    sub.x.rI2K <- function(x.=x, rI2K){
        if(length(dim(x))<3){
            return(x[rI2K[[1]], rI2K[[2]], drop=FALSE])
        } else {
            return(x[rI2K[[1]], rI2K[[2]],, drop=FALSE ])
        }
    }
    if(is.list(returnIndices2Keep)){
        if(length(returnIndices2Keep) != 2){
            stop('If returnIndices2Keep is a list, it must ',
                 'have length 2;  length = ',
                 length(returnIndices2Keep))
        }
        for(i in 1:2){
#           are returnIndices2Keep legal?
            ri <- returnIndices2Keep[[i]]
            oopsi <- which(abs(ri) > dim(x)[i])
            if(length(oopsi)>0){
                stop('returnIndices2Keep[[', i,
                     ']] has elements > dim(x)[', i, '].')
            }
            n0 <- sum(ri==0)
            if(n0>0){
                stop('returnIndices2Keep[[', i,
                     ']] includes 0s;  not allowed')
            }
            npos <- sum(ri>0)
            nneg <- sum(ri<0)
            if(npos*nneg>0){
                stop('returnIndices2Keep[[', i,
                     ']] includes both positive and ',
                     'negative numbers;  not allowed')
            }
            fp <- (ri%%1)
            if(any(fp>0)){
                stop('returnIndices2Keep[[', i,
                     ']] incudes non-integers; not allowed')
            }
        }
        return(sub.x.rI2K(rI2K=returnIndices2Keep))
    }
##
## 3.  !is.logical(returnIndices2Keep)?
##
    if(!is.logical(returnIndices2Keep)){
        stop('returnIndices2Keep must be either logical or a list;',
             '  class = ', class(returnIndices2Keep) )
    }
##
## 4.  compute indices2Keep
##
    maxAbs <- function(x){
        x.na <- is.na(x)
        if(na.rm){
            if(all(x.na)){
                return(-Inf)
            } else {
                return(max(abs(x[!x.na])))
            }
        } else {
            if(any(x.na)){
                return(Inf)
            } else {
                return(max(abs(x)))
            }
        }
    }
    max1 <- apply(x, 1, maxAbs)
    max2 <- apply(x, 2, maxAbs)
    max1p <- (max1>max2trim)
    max2p <- (max2>max2trim)
    imax1 <- which(max1p)
    imax2 <- which(max2p)
#
    if(length(imax1)>0){
        i1r <- range(imax1)
        i1 <- i1r[1]:i1r[2]
    } else i1 <- integer(0)
    if(length(imax2)>0){
        i2r <- range(imax2)
        i2 <- i2r[1]:i2r[2]
    } else i2 <- integer(0)
    indices2Keep <- list(index1=i1, index2=i2)
##
## 5.  return either indices or subsetted x
##
    if(returnIndices2Keep){
        return(indices2Keep)
    }
    return(sub.x.rI2K(rI2K=indices2Keep))
}
