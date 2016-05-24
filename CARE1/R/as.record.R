as.record <-
function(x){
t <- ncol(x) 
M <- nrow(x)
Mat <- matrix(0, ncol = t, nrow = 2 ^ t - 1)
for(i in 1:(2 ^ t - 1)) Mat[i,] <- rev(as.integer(intToBits(i))[1:t])
Lable <- capture.output(for(i in 1:(2^t-1)) cat(Mat[i,], "\n", sep=""))
z <- rep(0, 2 ^ t - 1)
for(i in 1 :(2 ^ t - 1))z[i] <- length(which(colSums(t(x) == Mat[i,]) == t))
names(z) <- Lable
return(z)
}
