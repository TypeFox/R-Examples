# partial derivatives of conditional correlation matrix w.r.t. its arguments
   vdR <- function(n){
     n1 <- n*(n-1)/2
     vecdR <- matrix(0, n1, n^2)
     M <- diag(-1,n)
     M[lower.tri(M)]<- 1:n1
     M <- M + t(M)
     for (i in 1:n1){
       m <- as.vector(M)
       m[m!=i] <- 0
       vecdR[i,] <- m
     }
     vecdR[vecdR!=0]<-1
     vecdR
   }
