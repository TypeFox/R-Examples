rinv <- function(dset){
nitem <- ncol(dset)-1
temp <- matrix(data = 0, nrow = nrow(dset), ncol = nitem, byrow = TRUE)
rorder <- matrix(data = 0, nrow = nrow(dset), ncol = nitem, byrow = TRUE)
rorder <- dset[,1:nitem]

## compute the ordering of the observations
for (i in 1:(nrow(rorder))){
for (j in 1:nitem){
for (k in 1:nitem){
if (rorder[i,j] == k){
temp[i,k] <- j
}
}
}
}

return(data.frame(temp,n=dset[,(nitem+1)]))
}