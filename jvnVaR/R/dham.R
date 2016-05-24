dham <-
function(x, y) {
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
dhamm[1] <- (-1 / v[1] + u / v[1] ^ 2) * v1[1]
dhamm[2] <- (-1 / v[1] + u / v[1] ^ 2) * v2[1]
dhamm[3] <- (-1 / v[1] + u / v[1] ^ 2) * v3[1]
for (i in 2:N) {
 v[i] <- x[3] + x[1] * u + x[2] * v[i-1]
v1[i] <- u + x[2] * v1[i-1]
v2[i] <- x[2] * v2[i-1] + v[i-1]
v3[i] <- 1 + x[2] * v3[i-1]
u <- y[i] ^ 2
dhamm[1] <- dhamm[1] + (-1 / v[i] + u / v[i] ^ 2) * v1[i]
dhamm[2] <- dhamm[2] + (-1 / v[i] + u / v[i] ^ 2) * v2[i]
dhamm[3] <- dhamm[3] + (-1 / v[i] + u / v[i] ^ 2) * v3[i]
}
return(dhamm)
}
