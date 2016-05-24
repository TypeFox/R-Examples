create.model <- function(observed,n){
#creates HWE model distribution from the matrix of observed genotype counts
obs <- as.matrix(observed,rownames.force=T)
r <- nrow(obs)

model <- mat.or.vec(r,r)
k <- 1
while(k<=r){
for(j in k:r){
if(j==k){
n_j <- 0
count <- 1
while(count <= j){
n_j <- n_j + obs[j,count]
count <- count + 1
}
count <- j
while(count <= r){
n_j <- n_j + obs[count,j]
count <- count + 1
}
model[j,k] <- (n_j)^2/(4*n)
}
if(j!=k){
n_j <- 0
count <- 1
while(count <= j){
n_j <- n_j + obs[j,count]
count <- count + 1
}
count <- j
while(count <= r){
n_j <- n_j + obs[count,j]
count <- count + 1
}
n_k <- 0
count <- 1
while(count <= k){
n_k <- n_k + obs[k,count]
count <- count + 1
}
count <- k
while(count <= r){
n_k <- n_k + obs[count,k]
count <- count + 1
}
model[j,k] <- ((n_j)*(n_k))/(2*n)
}
}
k <- k + 1
}

return(model)
}
