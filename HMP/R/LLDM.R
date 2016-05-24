LLDM <-
function(data, alpha){
N <- rowSums(data)
lprob <- vector(length=nrow(data))

for(i in 1:nrow(data)){
lprob[i] <- log(N[i]) + lbeta(sum(alpha), N[i]) -
sum(log(data[i,])) - sum(lbeta(alpha , as.numeric(data[i,])))
}

return(sum(lprob))
}
