##########################
####### solve3.QP ########
##########################

solve3.QP <- function(D, d, A, b, meq, factorized=FALSE) {
#### This code works for QP problem: min 1/2x'Dx-d'x
#### where the matrix D is diagonal and the constraints
#### are all equalities, i.e. t(A)x=b.
#### Inputs:
#### D should be a vector of length n, this means the matrix diag(D), but
#### if factorized=TRUE, D actually is diag(D)^(-1/2).
#### d is a vector of length n
#### A is a matrix of n x p
#### b is a vector with length p.  Finally meq = integer p.
#### The input meq are here for the compatibility with solve.QP in R package
#### quadprog. Written by Mai Zhou (mai@ms.uky.edu) Jan.30, 2001
D <- as.vector(D)
if(length(b)!=meq) stop("length of constraints not matched")
if(length(D)!=length(d)) stop("dimention of D and d not match")
if(dim(A)[1]!=length(D)) stop("dimention of D and A not match")
if(dim(A)[2]!=meq) stop("dimention of A not match with meq")

         if(!factorized) { D <- 1/sqrt(D) }
         QRout <- qr(D*A)
         temp <- rep(0,meq)
         if(any(b!=0)) {temp<-forwardsolve(t(qr.R(QRout)),b)}
         temp2 <- temp - t(qr.Q(QRout)) %*% (D * d)
         eta <- base::backsolve(qr.R(QRout), temp2)    #### added 2014/1/16 because SparseM also have a backsolve
         sol <- D^2 * (d + A %*% eta )
list( solution=sol )
}
