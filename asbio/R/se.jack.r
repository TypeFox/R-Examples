se.jack <- function(x){
n <- length(x)
x <- as.matrix(x)
sqrt(((n - 1)/n) * sum((x - mean(x))^2))
}

se.jack1 <- function(x){
n <- length(x)
x <- as.matrix(x)
sqrt(sum((x - mean(x))^2)/(n*(n - 1)))
}
