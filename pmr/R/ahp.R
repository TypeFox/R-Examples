ahp <- function(dset, sim_size=500){

nitem <- ncol(dset)

## check square
if (ncol(dset) != nrow(dset)){
message("Not a square matrix")
}

## check a_{ij} = 1/a_{ij}
recip <- 0

for (i in 1:nitem){
for (j in 1:nitem){
if (dset[i,j] < 1/dset[j,i]*0.99 || dset[i,j] > 1/dset[j,i]*1.01){
recip <- 1
}
}
}

if (recip == 1){
message("Please check a_{ij} = 1/a_{ij}")
}

## start if square matrix and a_{ij} = 1/a_{ij}
if (ncol(dset) == nrow(dset) && recip==0){

## Eigenvector Method
## weight <- Re(eigen(dset)$vector[1,])
## nth root of the product of the values are used to estimate the eigenvalue
weight <- 1:nitem
for (i in 1:nitem){
weight[i] <- prod(dset[i,])^(1/nitem)
}
temp_sum <- sum(weight)
for (i in 1:nitem){
weight[i] <- weight[i]/temp_sum
}
lambda_max <- Re(eigen(dset)$values[1])
CI <- (lambda_max-nitem)/(nitem-1)

## Saaty's inconsistency 
RI <- rep(1,sim_size)
for (i in 1:sim_size){
test <- matrix(data = 0, nrow = nitem, ncol = nitem, byrow = T)

## randomly drawn
weight_random <- sample(1:9,nitem*nitem,replace=T)
inverse_random <- sample(0:1,nitem*nitem,replace=T)

for (j in 0:(nitem*nitem-1)){
if ((j%/%nitem) < (j%%nitem)){
test[(j%/%nitem+1),(j%%nitem+1)] <- weight_random[j]
if (inverse_random[j] == 1){
test[(j%/%nitem+1),(j%%nitem+1)] <- 1/test[(j%/%nitem+1),(j%%nitem+1)]
}
}

if ((j%/%nitem) == (j%%nitem)){
test[(j%/%nitem+1),(j%%nitem+1)] <- 1
}

if ((j%/%nitem) > (j%%nitem)){
test[(j%/%nitem+1),(j%%nitem+1)] <- 1/test[(j%%nitem+1),(j%/%nitem+1)]
}

}

## store in RI
lambda_max_temp <- Re(eigen(test)$values[1])
CI_temp <- (lambda_max_temp-nitem)/(nitem-1)
RI[i] <- CI_temp
}

CR <- CI/mean(RI)

weight_RI <- Re(eigen(test)$vector[1,])
lambda_max_RI <- Re(eigen(test)$values[1])
RI[i] <- (lambda_max_RI-nitem)/(nitem-1)

## Koczkodaj's inconsistency

Ix <- 0

for (i in 1:(nitem-2)){
for (j in (i+1):(nitem-1)){
for (k in (i+2):nitem){
if (min(abs(1-dset[i,k]/dset[i,j]/dset[j,k]),abs(1-dset[i,j]*dset[j,k]/dset[i,k]))>Ix) {
Ix <- min(abs(1-dset[i,k]/dset[i,j]/dset[j,k]),abs(1-dset[i,j]*dset[j,k]/dset[i,k]))
}
}
}
}

## output
lst <- list(weighting=weight, Saaty=CR, Koczkodaj=Ix)
message("Summary of pairwise comparison matrics: ")
message("$weighting: weights of items; $Saaty: Saaty's inconsistency; $Koczkodaj: Koczkodaj's inconsistency")
return(lst)

}
}