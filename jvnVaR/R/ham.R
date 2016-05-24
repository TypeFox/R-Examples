ham <-
function(x, y) {
N <- length(y)
v <- c()
v[1] <- x[3] / (1 - x[1] - x[2])
u <- y[1] ^ 2
hamm <- -log(v[1]) - u / v[1]
for (i in 2:N) {
v[i] <- x[3] + x[1] * u + x[2] * v[i-1]
u <- y[i] ^ 2
hamm <- hamm - log(v[i]) - u / v[i]
}
return(hamm)
}
