
## Non-negative least squares with a trimmed W.
 
 ## (To be finished)

 ## @param H Left matrix H
 ## @param A Original matrix A
 ## @param zeta A logic matrix of the same dimension as A. TRUE -> cells to keep. FALSE -> cells to trim.
 ## @param alpha alpha
 ## @param k Reduced dimension.
 
 ## @return An updated matrix W.

 ## @examples
 ## #R code here showing how your function works (To be finished)

Nnls.trimW <- function(H, A, zeta, alpha, k, p)
  {
    ## Solution 2: if entire row is trimmed, randomly generate a row.
    fun1 = function(i, H, A, zeta){
        if(all(!c(zeta[i,]))){
            ##browser()
            warning(paste("Row",i,"is trimmed entirely in one of the iterations. It is temporarily replaced by random numbers from 0 to 1 to fit H."))
            x.sub = runif(ncol(A), 0, 1)
            H.ext = rbind(t(H), sqrt(alpha) * diag(k))
            x.ext = c(x.sub, rep(0,k))
            return(nnls(H.ext, x.ext)$x)
        }else{
            col.keep = c(zeta[i, ])
            H.trim = H[, col.keep]
            x.trim = A[i, col.keep]
            H.ext = rbind(t(H.trim), sqrt(alpha) * diag(k))
            x.ext = c(x.trim, rep(0,k))
            return(nnls(H.ext, x.ext)$x)
        }
    }
    
    ## ## Solution 3: if entire row is trimmed, randomly select 50% of the row to keep.
    ## fun1 = function(i, H, A, zeta){
    ##     if(all(!c(zeta[i,]))){
    ##         ##browser()
    ##         col.keep = c(zeta[i,])
    ##         col.keep[sample(1:length(col.keep), round(0.5 * length(col.keep)), replace = FALSE )] = TRUE
    ##         H.trim = H[,col.keep]
    ##         x.trim = A[i,col.keep]
    ##         H.ext = rbind(t(H.trim), sqrt(alpha) * diag(k))
    ##         x.ext = c(x.trim, rep(0,k))
    ##         return(nnls(H.ext, x.ext)$x)
    ##     }else{
    ##         ##browser()
    ##         col.keep = c(zeta[i,])
    ##         H.trim = H[,col.keep]
    ##         x.trim = A[i,col.keep]
    ##         H.ext = rbind(t(H.trim), sqrt(alpha) * diag(k))
    ##         x.ext = c(x.trim, rep(0,k))
    ##     return(nnls(H.ext, x.ext)$x)
    ##     }
    ## }
    return(t(sapply(1:p, fun1, H, A, zeta)))
}
