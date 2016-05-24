.packageName <- "sparseLDA"

normalizetest <- function(Xtst,Xn){
# normalize the columns of the test set Xtst using the mean and lengths of 
# the training set with output Xn from the function normalize.

ntst <- dim(Xtst)[1]
p <- dim(Xtst)[2]
m <- rep(Xn$mx,ntst)
dim(m) <- c(p,ntst)
Xtst <- Xtst-t(m)
v <- rep(Xn$vx,ntst)
dim(v) <- c(p,ntst)
v <- t(v)
Xtst <- Xtst[,Xn$Id]/v[1:ntst,Xn$Id]
}
 
