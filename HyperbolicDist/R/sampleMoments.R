### Functions to calculate the sample skewness and kurtosis
### Taken from the package e1071
skewness <- function (x, na.rm = FALSE){
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  sum((x - mean(x))^3)/(length(x)*sd(x)^3)
}

kurtosis <- function (x, na.rm = FALSE){
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  sum((x - mean(x))^4)/(length(x)*var(x)^2) - 3
}
