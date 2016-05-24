ATilde <- function(A){
 n <- NROW(A) 
 total.mean <- sum(A)/(n^2)
 m1 <- sapply(1:n, FUN=function(j) mean(A[j,]))
 m2 <- sapply(1:n, FUN=function(j) mean(A[,j]))
 A2 <- sapply(1:n,1:n,FUN=function(i,j) (A[i,j]-m1[i]-m2[j]+total.mean))
 return(A2)
}


