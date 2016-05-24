skewness <- function (x)
mean((x - mean(x))^3)/sd(x)^3

