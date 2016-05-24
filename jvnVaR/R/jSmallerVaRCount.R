jSmallerVaRCount <-
function(s,VaR){
n <- length(s)
v <- c()
j <- 1
for (i in 1:n) {
if (s[i] <= -VaR) {
v[j] <- i
j <- j + 1
}
}
k <- length(v)
for (i in k:2) {
v[i] <- v[i] - v[i-1]
}
return(v)
}
