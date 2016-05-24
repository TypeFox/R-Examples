estN.para <-
function(z, nhat){
t <- as.integer(log(length(z) + 1) / log(2))
Mat <- matrix(0, ncol = t, nrow = 2 ^ t - 1)
for(i in 1:(2 ^ t - 1)) Mat[i,] <- rev(as.integer(intToBits(i))[1:t])
ni <- sapply(1:t, function(k) Sub.n(z, t, Mat, k))
ui <- ni / nhat
names(ui) <- sapply(1:t, function(k) paste("u", k, sep=""))
rij <- matrix(0, t, t)
label <- matrix(0, t, t)

for(i in 1:(t - 1)){
for(j in (i + 1):t){
ni <- sum(z[which(Mat[,i] == 1)])
nj <- sum(z[which(Mat[,j] == 1)])
nij <- sum(z[which(Mat[,i] == 1 & Mat[,j] == 1)])
rij[i,j] = nhat * nij / (ni * nj) - 1
label[i,j] <- paste("r", i, j, sep="")
}
}
rij <- rij[upper.tri(rij)]
names(rij) <-  label[upper.tri(label)]
para <- c(ui, rij)
return(para)
}
