haml <-
function(x, y, L, U) {
N <- length(y)
v <- c()
v[1] <- x[3] / (1 - x[1] - x[2])
u <- y[1] ^ 2
if (y[1] <= L){
temp <- log(jNormCdf(L/sqrt(v[1])))
}
if (y[1] >= U){
temp <- log(jNormCdf(-U/sqrt(v[1])))
}
if (y[1] > L & y[1] < U){
temp <- - log(v[1]) - u / v[1]
}
hamm <- temp
for (i in 2:N) {
v[i] <- x[3] + x[1] * u + x[2] * v[i-1]
u <- y[i] ^ 2
if (y[i] <= L){
temp <- log(jNormCdf(L/sqrt(v[i])))
}
if (y[i] >= U){
temp <- log(jNormCdf(-U/sqrt(v[i])))
}
if (y[i] > L & y[i] < U){
temp <- - log(v[i]) - u / v[i]
}
hamm <- hamm + temp
}
return(hamm)
}
