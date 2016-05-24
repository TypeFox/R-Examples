library("RUnit")
library("kyotil")

test.matrix2 <- function() {

tolerance=1e-3
if(file.exists("D:/gDrive/3software/_checkReproducibility")) tolerance=1e-6
RNGkind("Mersenne-Twister", "Inversion")

d=1:4
X=matrix(1:12,4,3)
checkTrue(all(tXDX(X, d) == t(X) %*% diag(d) %*% X))

d1=1:3
d2=4:6
X=matrix(1:9,3,3)
checkTrue(all(DXD(d1, X, d2) == diag(d1) %*% X %*% diag(d2)))

S=matrix(c(1,2,3,2,4,5,3,5,8),3,3)
X=matrix(1:9,3,3)
checkTrue(all( symprod(S, X) == S %*% X))

x=1:3
y=4:6
S=matrix(c(1,2,3,2,4,5,3,5,8),3,3)
checkTrue(txSy(x, S, y) == drop(t(x)%*%S%*%y))


}


##check speed
#n=1e3
#d1=1:n
#d2=1:n
#X=matrix(1:n^2,n,n)
#all(DXD(d1, X, d2) == diag(d1) %*% X %*% diag(d2))
#system.time({DXD(d1, X, d2)})
##   user  system elapsed 
##   0.12    0.00    0.12 
#system.time({diag(d1) %*% X %*% diag(d2)})
##   user  system elapsed 
##   4.62    0.02    4.64 

#n=1e4
#x=1:n
#y=1:n+10
#S=diag(n)
#checkTrue(txSy(x, S, y) == drop(t(x)%*%S%*%y))
#system.time({txSy(x, S, y)})
##  user  system elapsed 
##   1.16    0.36    1.51 
#system.time({drop(t(x)%*%S%*%y)})
##   user  system elapsed 
##   1.95    0.00    1.95
