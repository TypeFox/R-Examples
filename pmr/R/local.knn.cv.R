local.knn.cv <- function(dset,covariate.test,covariate,cv=10,k.max=20,method.cv="mean"){

nitem <- ncol(dset)
ncov <- ncol(covariate)

## error function = Kendall's distance

rdist <- function(x,y){
d <- 0
for (j in 1:(nitem-1)){
for (k in (j+1):nitem){
if ((x[j]-x[k])*(y[j]-y[k]) < 0) {d <- d+1}
}
}
d
}

## create cv dataset for cross validation

cv_label <- sample(1:nrow(dset),nrow(dset),replace=F)
cv_label <-cv_label%%cv+1

## test k

error_k <- rep(0,k.max)

for (i in 1:k.max){

## test cv

error_cv <- rep(0,cv)
for (j in 1:cv){

## extract training ranking, training covariate, and testing covariate
dset_cv <- matrix(data = 0, nrow = nrow(dset)-table(cv_label)[j], ncol = nitem, byrow = TRUE)
covariate_cv <- matrix(data = 0, nrow = nrow(dset)-table(cv_label)[j], ncol = ncov, byrow = TRUE)
covariate.test_cv <- matrix(data = 0, nrow = table(cv_label)[j], ncol = ncov, byrow = TRUE)
dset.test_cv <- matrix(data = 0, nrow = table(cv_label)[j], ncol = nitem, byrow = TRUE)

count_cv1 <- 1
count_cv2 <- 1
for (k in 1:nrow(dset)){
if (cv_label[k] != j){
for (m in 1:nitem){
dset_cv[count_cv1,m] <- dset[k,m]
covariate_cv[count_cv1,m] <- covariate[k,m]
}
count_cv1 <- count_cv1 + 1
}	
else{
for (m in 1:nitem){	
dset.test_cv[count_cv2,m] <- dset[k,m]
covariate.test_cv[count_cv2,m] <- covariate[k,m]
}
count_cv2 <- count_cv2 + 1	
}
}

## fit knn to dset_cv

dset.test_knn <- local.knn(dset_cv,covariate.test_cv,covariate_cv,knn.k=i,method=method.cv)

## compute error
for (k in 1:nrow(dset.test_cv)){
error_cv[j] <- error_cv[j] + rdist(dset.test_cv[k,],dset.test_knn[k,])
}

} ## fin test cv

## compute error for this k
error_k[i] <- mean(error_cv)
} ## fin test k

message("k=",which.min(error_k))
return(local.knn(dset_cv,covariate.test_cv,covariate_cv,knn.k=which.min(error_k),method=method.cv))

}