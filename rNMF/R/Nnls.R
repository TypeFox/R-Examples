## Performs Non-negative least square regression on a matrix.
## 
## Given two matrices M, C >= 0, find W >= 0 such that ||C - MW||_2^2 is minimized.
## (This parts to be finished)
##
## @param M Design matrix. \code{M}
## @param C Response matrix \code{C}
## 
## @return output W, a non-negative matrix (To be finished)
##
## @keywords keywords
##
## @export
## 
## @examples
## #R code here showing how your function works (To be finished)

Nnls = function(M = M, C = C){
    fun1 = function(x){ return(nnls(M, x)$x)}
    return(apply(C, 2, fun1))
}
