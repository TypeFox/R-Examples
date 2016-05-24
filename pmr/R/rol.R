rol <- function(dset,covariate){

nitem <- ncol(dset)
dset <- cbind(dset,rep(1,nrow(dset)))
oset <- rinv(dset[,1:(nitem+1)])
oset <- oset[,1:ncol(oset)-1]
covariate <- cbind(rep(1,nrow(oset)),covariate)
ncov <- ncol(covariate)

names(covariate)[1] <- "Intercept"

## loglikelihood function
loglik_rol <- function(slope1){
## transform linear parameter into matrix
slope <- matrix(data = 0, nrow = ncov, ncol = nitem-1, byrow = TRUE)
for (i in 1:ncov){
for (j in 1:nitem-1){
slope[i,j] <- slope1[(i-1)*(nitem-1)+j]
}
}

util <- matrix(data = 0, nrow = nrow(oset), ncol = nitem-1, byrow = TRUE)
## reference item: last item
util_last <- rep(0,nrow(oset))

for (i in 1:nrow(oset)){
for (j in 1:nitem-1){
for (k in 1:ncov){
util[i,j] <- util[i,j] + covariate[i,k] * slope[k,j]
}
}
}

util_combine <- cbind(util, util_last)

pr <- rep(1,nrow(oset))

for (i in 1:nrow(oset)){
deno <- sum(exp(util_combine[i,]))
for (j in 1:(nitem-1)){
pr[i] <- pr[i] * exp(util_combine[i,oset[i,j]]) / deno
deno <- deno - exp(util_combine[i,oset[i,j]])
}
}

## ll <-sum(-log(pr))

ll <- 0
for (i in 1:nrow(oset)){
if (pr[i]>0.0000000000000001 & !is.na(pr[i])){
ll <- ll - log(pr[i])
}
else{}
}

ll
}

## out1 <- optim(rep(0,ncov*nitem-ncov), loglik_rol, NULL, method = "BFGS", hessian = TRUE)
out1 <- nlm(loglik_rol, rep(0,ncov*nitem-ncov), hessian=TRUE)
##require(stats4)
message("Maximum Likelihood Estimation of the Rank-ordered Logit Model")
out2 <- new("mle")
out2@coef <- out1$estimate
out2@fullcoef <- out1$estimate
out2@vcov <- solve(out1$hessian)
out2@min <- out1$minimum
out2@details <- out1
out2@minuslogl <- loglik_rol
out2@method <- "Newton-Raphson"

## labeling of coef
for (i in 1:ncov){
for (j in 1:nitem-1){
names(out2@coef)[(i-1)*(nitem-1)+j] <- paste("Beta",i-1,"item",j, sep = "")
names(out2@fullcoef)[(i-1)*(nitem-1)+j] <- paste("Beta",i-1,"item",j, sep = "")
}
}

return(out2)
}