.spur<-function(U){
    sum(diag(U))
}

.orthComplement<-function(W) {
    ##orthogonal complement of <W>: <W>orth=<Worth>
    rW    <- rankMatrix(W)
    Worth <- qr.Q(qr(cbind(W)),
                  complete=TRUE)[,-c(1:rW),drop=FALSE]
    Worth
}





## Old UHH-code below
## Completely obsolete

## .spurAB<-function(A,B){
##   sum(A*t.default(B))
## }
## # if A eller B er symmetrisk så er trace(AB)=sum(A*B)


## .matrixNullSpace<-function(B,L) {
##   ## find A such that <A>={Bb| b in Lb=0}
##   if ( ncol(B) != ncol(L) ) {
##      stop('Number of columns of B and L not equal \n')
##   }
##   A <- B %*% .orthComplement(t(L))
##   # makes columns of  A orthonormal:
##   A <- qr.Q(qr(A))[,1:rankMatrix(A)]
##   A
## }


## .colSpaceCompare<-function(X1,X2) {
##   ## X1 X2: matrices with the ssme number of rows
##   ## results r  (Ci column space of Xi)
##   ## r=1 C1 < C2
##   ## r=2 C2 < C1
##   ## r=3 C1==C2
##   ## r=4 C1 intersect C2 NOTempty but neither the one contained in the other
##   ## r=5 C1 intersect C2 = empty
##   if (nrow(X1)!= nrow(X2)){
##     stop("\n number of rows of X1 and X2 must be equal") }
##   r1 <-rankMatrix(X1)
##   r2 <-rankMatrix(X2)
##   r12<-rankMatrix(cbind(X1,X2))
##   r  <-
##     if (r12 <= pmax(r1,r2)) {
##       if (r1==r2) 3 else {
##         if (r1>r2) 2 else 1
##       }
##     } else {
##       if (r12==(r1+r2))  5 else 4
##     }
##   r
## }
