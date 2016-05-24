
## Non-negative least squares with a trimmed H.
 
## @param W: Left matrix 
## @param A: Original matrix A
## @param zeta: A logic matrix of the same dimension as A. TRUE -> cells to keep. FALSE -> cells to trim.
## @param alpha alpha
## @param k Reduced dimension.
## @param n The number of columns of A.
 
## @return An updated matrix W.

## @examples
## #R code here showing how your function works (To be finished)

## Nnls.trimH for variation "cell"

Nnls.trimH = function(W, A, zeta, beta, k, n){
    fun1 = function(j, W, A, zeta){
        if(all(!c(zeta[,j]))){
            stop("OK, an entire column is trimmed. I really don't know what to do.")
        }else{
            row.keep = c(zeta[,j])
            W.trim = W[row.keep,]
            x.trim = A[row.keep,j]
            W.ext = rbind(W.trim, sqrt(beta) * matrix(1,1,k))
            x.ext = c(x.trim, 0)
        }
        return(nnls(W.ext,x.ext)$x)
    }
    return(sapply(1:n, fun1, W, A, zeta))
}
