rankagg <- function(dset){

nitem <- ncol(dset)

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
if (temp_ind == 0) {n[j] <- n[j] + 1}
}
}
test2 <- cbind(test, n)
rankagg <- test2

# reduced those with frequency = 0

count_nonzero <- 0
for (i in 1:factorial(nitem)){
if (test2[i,nitem+1] != 0) {count_nonzero <- count_nonzero + 1}
}
# created a reduced dataset
test3 <- matrix(data = 0, nrow = count_nonzero, ncol = nitem, byrow = TRUE)
n <- rep(0,count_nonzero)
temp_ind <- 1
for (i in 1:factorial(nitem)){
if (test2[i,nitem+1] != 0) {
n[temp_ind] <- n[temp_ind] + test2[i,nitem+1]
for (j in 1:nitem){
test3[temp_ind,j] <- test2[i,j]	
}
temp_ind <- temp_ind + 1
}
}
rankagg <- cbind(test3,n)	


}