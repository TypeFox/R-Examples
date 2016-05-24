## ## This function reorders (permutes) two matrices A and B which have the same shape so that their column structures are as similar as possible. The output value is the permutation order for the second matrix.
## ## WARNING: this function is slow!!
## ## NOTE: Currently not enabled in the main function "rNMF".
## per = function(A,B){
##     fun1 = function(x){
##         sum(apply(A - B[,x], 2, l2))
##     }
##     m = ncol(A)
##     per = permn(m)
##     diff = lapply(per, fun1)
##     return(per[[which.min(diff)]])
## }

## "find.row" takes the (x,y) coordinates in A in the format xypair = c(x,y), and returns the corresponding row index of v.
## "find.x" takes the index of v, and returns the x coordinator in A.
## "find.y" takes the index of v, and returns the y coordinator in A.
find.row = function(xypair, p){
    xypair = unlist(xypair)
    return(p * (xypair[2] - 1) + xypair[1])
}
