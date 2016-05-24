simule_MC <-
function(transmat,prior,T) {

y = matrix(0,1,T)
y[1] = which(rmultinom(1, size = 1, prob = prior)==1)
for(i in 2:T) {
y[i]=which(rmultinom(1, size = 1, prob = transmat[y[i-1], ])==1)
}
return(y)
}
