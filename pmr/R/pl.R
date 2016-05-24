pl <- function(dset){

oset <- rinv(dset)
nitem <- ncol(dset)-1

## loglikelihood function
loglik_pl <- function(util){
pr <- rep(1,nrow(oset))

temp_mean <- mean(util)
for (i in 1:nitem){
util[i] <- util[i] - temp_mean
}

for (i in 1:nrow(oset)){
deno <- sum(exp(util))
for (j in 1:(nitem-1)){
pr[i] <- pr[i] * exp(util[oset[i,j]]) / deno
deno <- deno - exp(util[oset[i,j]])
}
}

ll <- rep(0,nrow(oset))
for (i in 1:nrow(oset)){
ll[i] = -log(pr[i])*oset[i,(nitem+1)]
}
sum(ll)
}

out1 <- optim(rep(1,nitem), loglik_pl, NULL, method = "BFGS", hessian = TRUE)
##out1 <- mle(minuslogl = loglik_pl, start = list(util=rep(1,nitem)), method = "BFGS", nobs=sum(dset[,nitem+1]))
##require(stats4)
test <- matrix(data = 0, nrow = factorial(nitem), ncol = nitem, byrow = TRUE)
temp1 <- 1:nitem
i <- 1

## generate a list of all possible rankings
for (j in 1:(nitem^nitem-1)){
temp1[1] <- nitem - j%%nitem
temp2 <- j - j%%nitem
for (k in nitem:2){
temp1[k] <- nitem - temp2%/%(nitem^(k-1))
temp2 <- temp2 - (nitem-temp1[k])*(nitem^(k-1))
}
temp2 <- 0
for (l in 1:nitem){
for (m in 1:nitem){
if (temp1[l] == temp1[m] && l != m){
temp2 <- 1
}
}
}
if (temp2 == 0){
for (p in 1:nitem){
test[i,p] = temp1[p]
}
i <- i+1
}
}

n <- rep(0,factorial(nitem))
for (j in 1:factorial(nitem)){
for (k in 1:nrow(dset)){
temp_ind <- 0
for (l in 1:nitem){
if (test[j,l] != dset[k,l]) {temp_ind <- temp_ind + 1}
}
if (temp_ind == 0) {n[j] <- dset[k,nitem+1]}
}
}
test2 <- cbind(test, n)
test2 <- rinv(test2)

## compute expected value
pro <- rep(1,factorial(nitem))
fitted <- rep(0,factorial(nitem))
for (i in 1:factorial(nitem)){
deno <- sum(exp(out1$par))
for (j in 1:(nitem-1)){
pro[i] <- pro[i] * exp(out1$par[test2[i,j]]) / deno
deno <- deno - exp(out1$par[test2[i,j]])
}
fitted[i] <- pro[i]*sum(n)
}

ss <- rep(0,factorial(nitem))
for (j in 1:factorial(nitem)){
ss[j] <- (fitted[j] - n[j])^2/fitted[j]
}

#lst <- list(loglik=out1$value, par=out1$par, se=(diag(solve(out1$hessian)))^0.5, fit.value=fitted, residual=sum(ss))
#return(lst)
#message("Modal ranking: ", modal)
message("Maximum Likelihood Estimation of the Luce Model")
message("Chi-square residual statistic: ", round(sum(ss), digits = 2), ", df: ", factorial(nitem))
out2 <- new("mle")
out2@coef <- out1$par
out2@fullcoef <- out1$par
out2@vcov <- solve(out1$hessian)
out2@min <- out1$value
out2@details <- out1
out2@minuslogl <- loglik_pl
out2@method <- "BFGS"
return(out2)
}
