local.knn <- function(dset,covariate.test,covariate,knn.k=1,method="mean"){

nitem <- ncol(dset)
ncov <- ncol(covariate)

## compute distance

## standardized covariates

covariate_z <- matrix(data = 0, nrow = nrow(dset), ncol = ncov, byrow = TRUE)
covariate.test_z <- matrix(data = 0, nrow = nrow(covariate.test), ncol = ncov, byrow = TRUE)
for (i in 1:ncov){
covariate_z[,i] <- (covariate[,i]-mean(covariate[,i]))/sd(covariate[,i])
covariate.test_z[,i] <- (covariate.test[,i]-mean(covariate[,i]))/sd(covariate[,i])
}

dist <- matrix(data = 0, nrow = nrow(covariate.test), ncol = nrow(dset), byrow = TRUE)
## compute distance for all judges
for (i in 1:nrow(covariate.test_z)){
for (j in 1:nrow(dset)){
for (k in 1:ncov){
dist[i,j] <- dist[i,j] + (covariate.test_z[i,k]-covariate_z[j,k])*(covariate.test_z[i,k]-covariate_z[j,k])
}
}
}

## max dist for i=j
##for (i in 1:nrow(dset)){
##for (j in 1:nrow(dset)){
##if (i == j){
##dist[i,i] <- max(dist[i,]) + 1 
##}
##}
##}

## compute ranking of distance for all observations

dist_rank <- matrix(data = 0, nrow = nrow(covariate.test), ncol = nrow(dset), byrow = TRUE)
for (i in 1:nrow(covariate.test)){
dist_rank[i,] <- rank(dist[i,],ties.method="random")
}

## create ranking dataset from kNN
dset_knn <- array(0, dim=c(nrow(covariate.test),knn.k,nitem))
for (i in 1:nrow(covariate.test)){
for (j in 1:nrow(dset)){
if (dist_rank[i,j] <= knn.k){
for (k in 1:nitem){
dset_knn[i,dist_rank[i,j],k] <- dset[j,k]
}
}
}
}

dset_predict <- matrix(data = 0, nrow = nrow(covariate.test), ncol = nitem, byrow = TRUE)

## mean rank
if (method=="mean"){
mean_knn <- matrix(data = 0, nrow = nrow(covariate.test), ncol = nitem, byrow = TRUE)
for (i in 1:nrow(covariate.test)){
for (j in 1:nitem){
mean_knn[i,j] <- mean(dset_knn[i,,j])
}
}

## assign rank, if tie mean rank then random
for (i in 1:nrow(covariate.test)){
dset_predict[i,] <- rank(mean_knn[i,],ties.method="random")
}

}

## luce model
if (method=="pl"){
oset_knn <- array(0, dim=c(nrow(covariate.test),knn.k,nitem))
temp <- matrix(data = 0, nrow = knn.k, ncol = nitem, byrow = TRUE)

for (i in 1:nrow(covariate.test)){
rorder <- matrix(data = 0, nrow = knn.k, ncol = nitem, byrow = TRUE)
rorder <- dset_knn[i,,]

## compute the ordering of the observations
for (h in 1:(nrow(rorder))){
for (j in 1:nitem){
for (k in 1:nitem){
if (rorder[h,j] == k){
temp[h,k] <- j
}
}
}
}

## fit Luce model on the ordering temp[]
loglik_pl <- function(util){
pr <- rep(1,nrow(temp))

temp_mean <- mean(util)
for (i in 1:nitem){
util[i] <- util[i] - temp_mean
}

for (i in 1:nrow(temp)){
deno <- sum(exp(util))
for (j in 1:(nitem-1)){
pr[i] <- pr[i] * exp(util[temp[i,j]]) / deno
deno <- deno - exp(util[temp[i,j]])
}
}

ll <- rep(0,nrow(temp))
for (i in 1:nrow(temp)){
ll[i] = -log(pr[i])
}
sum(ll)
}

out1 <- optim(rep(1,nitem), loglik_pl, NULL, method = "BFGS", hessian = TRUE)

## assign rank by Luce parameter, if tie mean rank then random
dset_predict[i,] <- rank(out1$par,ties.method="random")
dset_predict[i,] <- nitem+1-dset_predict[i,]

}
}

return(dset_predict)

}