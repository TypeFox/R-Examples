dhaml <-
function(x, y, L, U) {
 N <- length(y)
 v <- c()
v1 <- c()
v2 <- c()
v3 <- c()
dhamm <- c()
 v[1] <- x[3] / (1 - x[1] - x[2])
v1[1] <- x[3] / (1 - x[1] - x[2]) ^ 2
v2[1] <- x[3] / (1 - x[1] - x[2]) ^ 2
v3[1] <- 1 / (1 - x[1] - x[2])
u <- y[1] ^ 2
if (y[1] <= L){
temp <- -L*jNormPdf(L/sqrt(v[1]))/(jNormCdf(L/sqrt(v[1]))*sqrt(v[1])^3)
}
if (y[1] >= U){
temp <- U*jNormPdf(-U/sqrt(v[1]))/(jNormCdf(-U/sqrt(v[1]))*sqrt(v[1])^3)
}
if (y[1] > L & y[1] < U){
temp <- -1 / v[1] + u / v[1] ^ 2
}
dhamm[1] <- (-1 / v[1] + u / v[1] ^ 2) * v1[1]
dhamm[2] <- (-1 / v[1] + u / v[1] ^ 2) * v2[1]
dhamm[3] <- (-1 / v[1] + u / v[1] ^ 2) * v3[1]
for (i in 2:N) {
 v[i] <- x[3] + x[1] * u + x[2] * v[i-1]
v1[i] <- u + x[2] * v1[i-1]
v2[i] <- x[2] * v2[i-1] + v[i-1]
v3[i] <- 1 + x[2] * v3[i-1]
u <- y[i] ^ 2
if (y[i] <= L){
temp <- -L*jNormPdf(L/sqrt(v[i]))/(jNormCdf(L/sqrt(v[i]))*sqrt(v[i])^3)
}
if (y[i] >= U){
temp <- U*jNormPdf(-U/sqrt(v[i]))/(jNormCdf(-U/sqrt(v[i]))*sqrt(v[i])^3)
}
if (y[i] > L & y[i] < U){
temp <- -1 / v[i] + u / v[i] ^ 2
}
dhamm[1] <- dhamm[1] + temp * v1[i]
dhamm[2] <- dhamm[2] + temp * v2[i]
dhamm[3] <- dhamm[3] + temp * v3[i]
}
return(dhamm)
}
