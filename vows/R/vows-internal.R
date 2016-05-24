#' Internal functions for the vows package
#' 
#' These functions are ordinarily not to be called by the user, but if you
#' contact the authors with any questions about them, we'll do our best to
#' clarify matters.
#' 
#' 
#' @aliases get.ind rlr.xz vec2img
#' @author Lei Huang, Yin-Hsiu Chen, Philip Reiss \email{phil.reiss@@nyumc.org} and Ruixin Tan
#' @name vows-internal
#' @keywords internal
#' @export get.ind 
#' @export rlr.xz
#' @export vec2img
NULL

#' @rdname vows-internal
get.ind <-
function(arr3d) {
    flag = rep(NA, dim(arr3d)[1])
    for (i in  1:length(flag))    {
        cat("First-dimension slice",i,"\n")
        flag[i] = any(arr3d[i,,]!=0)
    }
    x.ind = which(flag)
    
    flag = rep(NA, dim(arr3d)[2])
    for (i in  1:length(flag))    {
        cat("Second-dimension slice",i,"\n")
        flag[i] = any(arr3d[,i,]!=0)
    }
    y.ind = which(flag)
	
    flag = rep(NA, dim(arr3d)[3])
    for (i in  1:length(flag))    {
        cat("Third-dimension slice",i,"\n")
        flag[i] = any(arr3d[,,i]!=0)
    }
    z.ind = which(flag)
    return(list(x.ind, y.ind, z.ind))
}

#' @rdname vows-internal
rlr.xz <-
function(x, nbasis=15, norder=4, nulldim=NULL, B=NULL, P) {
    if (is.null(B))   {
        bsb = create.bspline.basis(range(x, na.rm=TRUE), nbasis, norder)
        B = eval.basis(x, bsb)
        if (is.null(nulldim)) nulldim = norder-2
        P = getbasispenalty(bsb, nulldim)
    }
    else nbasis = ncol(B)
    svd.pen = svd(P)
    Us = svd.pen$u[ , 1:(nbasis-nulldim)]; Un = svd.pen$u[ , (nbasis-nulldim+1):nbasis]  
    X = B %*% Un
    Z = scale(B %*% Us, FALSE, sqrt(svd.pen$d[1:(nbasis-nulldim)]))  # see Wood (2004)
    list(X=X, Z=Z)
}

#' @rdname vows-internal
vec2img<-function(vec,logicarray){
    if(length(vec)!=sum(logicarray)) stop("dimension not fit")
    trueinds<-which(logicarray==TRUE,arr.ind=TRUE)
    result<-array(NA,dim(logicarray))
    result[logicarray]<-vec
    return(result)
}

