# prepare random sparse matrix
library(Matrix)
library(rmumps)
n=2000
a=Matrix(0, n, n)
set.seed(7)
ij=sample(1:(n*n), 15*n)
a[ij]=runif(ij)
diag(a)=0
diag(a)=-rowSums(a)
a[1,1]=a[1,1]-1
am=Rmumps$new(a)
b=as.double(a%*%(1:n)) # rhs for an exact solution vector 1:n
# following time includes symbolic analysis, LU factorization and system solving
system.time(x<-am$solve(b))
bb=2*b
# this second time should be much shorter
# as symbolic analysis and LU factorization are already done
system.time(xx<-am$solve(bb))
# compare to Matrix corresponding times
system.time(xm<-solve(a, b))
system.time(xxm<-solve(a, bb))
# compare to Matrix precision
range(x-1:n)  # mumps
range(xm-1:n) # Matrix

# matrix inversion
system.time(aminv <- am$inv())
system.time(ainv <- solve(a)) # the same in Matrix

# clean up by hand to avoid possible interference between gc() and
# Rcpp object destructor after unloading this namespace
rm(am)
gc()
