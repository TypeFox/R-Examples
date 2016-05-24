ddhame <-
function(x, y) {
 N <- length(y)
 v <- c()
v1 <- c()
v11 <- c()
 v[1] <- y[1] ^ 2
v1[1] <- 0
u <- y[1] ^ 2
v11[1] <- 0
ddhamm <- (-1 / v[1] + u / v[1] ^ 2) * v11[1] + (1 / v[1] ^ 2 - 2 * u / v[1] ^ 3) * v1[1] ^ 2
for (i in 2:N) {
  v[i] <- (1-x) * u + x * v[i-1]
 v1[i] <- -u + x * v1[i-1]
v11[i] <- x * v11[i-1] + v1[i-1]
     u <- y[i] ^ 2
ddhamm <- ddhamm + (-1 / v[i] + u / v[i] ^ 2) * v11[i] + (1 / v[i] ^ 2 - 2 * u / v[i] ^ 3) * v1[i] ^ 2
}
return(ddhamm)
}
