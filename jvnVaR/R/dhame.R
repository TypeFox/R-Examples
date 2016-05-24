dhame <-
function(x, y) {
 N <- length(y)
 v <- c()
v1 <- c()
 v[1] <- y[1] ^ 2
v1[1] <- 0
u <- y[1] ^ 2
dhamm <- (-1 / v[1] + u / v[1] ^ 2) * v1[1]
for (i in 2:N) {
 v[i] <- (1-x) * u + x * v[i-1]
v1[i] <- -u + x * v1[i-1]
u <- y[i] ^ 2
dhamm <- dhamm + (-1 / v[i] + u / v[i] ^ 2) * v1[i]
}
return(dhamm)
}
