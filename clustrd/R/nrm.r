nrm <- function(A){
# makes columnwise normalized version of A
n = dim(A)[1]
m = dim(A)[2]

d = apply(A^2,2,sum)^.5
N = A / (as.matrix(rep(1,n))%*%d)
}
