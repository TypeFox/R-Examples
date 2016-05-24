estN.HSC <-
function(z){
t <- as.integer(log(length(z) + 1) / log(2))
Mat <- matrix(0, ncol = t, nrow = 2 ^ t - 1)
for(i in 1:(2 ^ t - 1)) Mat[i,] <- rev(as.integer(intToBits(i))[1:t])
M <- sum(z)
D <- M - mean(sapply(1:t, function(j) Sub.S(z, t, Mat, j)))
Chat <- 1 - mean(sapply(1:t, function(j) Sub.S(z, t, Mat, j)) / sapply(1:t, function(j) Sub.n(z, t, Mat, j)))

MatA <- matrix(0, t, t)
MatB <- matrix(0, t, t)
MatC <- matrix(0, t, t)

for(i in 1:(t - 1)){
ni <- Sub.n(z, t, Mat, i)
for(j in (i + 1):t){
nj <- Sub.n(z, t, Mat, j)
MatA[i, j] = Sub.A(z, t, Mat, i, j)
MatB[i, j] = Sub.B(z, t, Mat, i, j)
MatC[i, j] = MatA[i, j] * MatB[i, j] /(ni * nj)
}
}

x <- t * D - sum(MatA)
y <- (t * D - sum(MatA)) / (t * Chat - sum(MatC))
Nhat=ifelse(x < 0.0001 || y <M, M, y)
names(Nhat) <- "Nhat"
return(Nhat)
}
