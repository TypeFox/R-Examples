mdpref <- function(dset,rank.vector=FALSE,ndim=2){

nitem <- ncol(dset)-1
nobs <- sum(dset[,nitem+1])
X <- matrix(data = 0, nrow = nobs, ncol = nitem, byrow = TRUE)

## check if dimension is smaller than or equal to the number of items
if (ndim > nitem){
message("Number of dimension has to be smaller than or equal to the number of items.")
}

if (ndim < nitem-1){
# enlarge ranking dataset
temp <- 1
for (i in 1:nrow(dset)){
if (dset[i,nitem+1] > 0){
for (j in 1:dset[i,nitem+1]){
for (k in 1:nitem){
X[temp,k] <- dset[i,k]
}
temp <- temp + 1
}
}
}

# reverse coding
for (i in 1:nrow(X)){
for (j in 1:nitem){
X[i,j] <- nitem-X[i,j]+1
}
}

# standardize
temp_mean <- rep(mean(X[1,]),nobs)
temp_sd <- rep(sd(X[1,]),nobs)

for (i in 1:nrow(X)){
for (j in 1:nitem){
X[i,j] <- (-X[i,j] + temp_mean[i])/temp_sd[i]
}
}

# multidimensional preference
md <- svd(X, nu=ndim, nv=ndim)
D <- matrix(data = 0, nrow = ndim, ncol = ndim, byrow = TRUE)
for (i in 1:ndim){
D[i,i] <- md$d[i]
}

# plot item graph
A_full <- md$u*(nitem-1)^0.5
B_full <- md$v %*% D/(nitem-1)^0.5
A_2 <- md$u[,1:2]*(nitem-1)^0.5
B_2 <- md$v[,1:2] %*% D[1:2,1:2]/(nitem-1)^0.5
label <- 1:nitem
plot(B_2, col="blue", xlab="Dimension 1", ylab="Dimension 2", type="n")
text(B_2[,1],B_2[,2],labels=label)

# draw rank vector
A_2 <- -A_2
A_full <- -A_full

d_info <- matrix(data = 0, nrow = ndim, ncol = ndim, byrow = TRUE)
for (i in 1:ndim){
for (j in 1:ndim){
d_info[i,j] <- max(B_full[,i])/max(A_full[,j])
}
}
mind <- min(d_info[1:2,1:2])
A_2_1 <- A_2 * mind
if (rank.vector==TRUE){
for (i in 1:nrow(A_2_1)){
testx <- c(0,A_2_1[i,1])
testy <- c(0,A_2_1[i,2])
lines(testx, testy)
}
}

coord <- matrix(data = 0, nrow = nrow(dset), ncol = nitem+1+ndim, byrow = TRUE)
# write output ranking
for (i in 1:nrow(dset)){
for (j in 1:(nitem+1)){
coord[i,j] <- dset[i,j]
}
}

temp <- 2
while (coord[temp,(nitem+1)] == 0){
temp <- temp + 1
}

coord[1,(nitem+2)] <- A_2_1[1,1]
coord[1,(nitem+3)] <- A_2_1[1,2]
for (i in 1:ndim){
coord[1,(nitem+1+i)] <- A_full[1,i]
}

while(temp < nrow(dset)){
for (i in 2:nrow(A_full)){
if (round(A_full[(i-1),1],digits=5) != round(A_full[i,1],digits=5) && round(A_full[(i-1),2],digits=5) != round(A_full[i,2],digits=5)){
for (j in 1:ndim){
coord[temp,(nitem+1+j)] <- A_full[i,j]
}
temp <- temp + 1
if (temp < nrow(dset)){
while (coord[temp,(nitem+1)] == 0){
temp <- temp + 1
}
}
}
}
}

var_explain_d <- 0
for (i in 1:ndim){
var_explain_d <- var_explain_d + D[i,i]
}

lst <- list(item=B_full, ranking=coord, explain=var_explain_d/sum(md$d))
return(lst)
}

}