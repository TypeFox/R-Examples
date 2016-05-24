#' Voxelwise linear models
#' 
#' This is a wrapper function for \code{\link{lm.mp}} to handle 3D image
#' responses.
#' 
#' 
#' @param arr4d a 4-dimensional response array, where the first 3 dimensions
#' refer to spatial coordinates and the last dimension corresponds to different
#' images.
#' @param formula,store.fitted see \code{\link{lm.mp}}.
#' @return An object of class \code{"\link{lm.mp}"}, with two changes. (1) If
#' \code{store.fitted = TRUE}, the fitted values are given as a 4-dimensional
#' array. (2) A \code{call} component is included.
#' @author Lei Huang \email{huangracer@@gmail.com}, Yin-Hsiu Chen
#' \email{enjoychen0701@@gmail.com}, and Philip Reiss
#' \email{phil.reiss@@nyumc.org}
#' @seealso \code{\link{lm.mp}}
#' @examples
#' 
#' data(test)
#' d4 = test$d4
#' x = test$x
#' lmobj = lm4d(d4, ~x)
#' 
#' # Convert d4 to a matrix, and confirm that lm.mp() gives the same results as lm4d()
#' d4.2 = d4
#' dim(d4.2) = c(prod(dim(d4)[1:3]), dim(d4)[4])
#' Y = t(d4.2)
#' lmobj2 = lm.mp(Y, ~x)
#' all.equal(lmobj$coef, lmobj2$coef)
#' @export
lm4d <- function(arr4d, formula, store.fitted=FALSE) {
    dim.4d = dim(arr4d)
    has.data<-attributes(arr4d)$has.data
    N=NROW(model.matrix(formula))
    dim(arr4d)=c(prod(dim(arr4d)[1:3]),N)
    Y.d = t(arr4d[(which(has.data)),])  
    rm(arr4d)  # i.e., remove from current environment
    Y.fit=lm.mp(Y=Y.d,formula=formula, store.fitted=store.fitted)
    if (store.fitted) {
        fit.value=array(NA,dim=dim.4d)
        for (j in 1:N) fit.value[,,,j][which(has.data, arr.ind=TRUE)]=Y.fit$fitted[j,]
        Y.fit$fitted=fit.value
    }
    Y.fit$call = match.call()
    return(Y.fit)
}
