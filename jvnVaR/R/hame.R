hame <-
function(x, y) {
N <- length(y)
v <- c()
v[1] <- y[1] ^ 2
u <- y[1] ^ 2
hamm <- -log(v[1]) - u / v[1]
for (i in 2:N) {
v[i] <- (1-x) * u + x * v[i-1]
u <- y[i] ^ 2
hamm <- hamm - log(v[i]) - u / v[i]
}
return(hamm)
}
