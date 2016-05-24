print.TOS2D <-function (x, ...) {

summary(x)

cat("Number of bootstrap realizations:",length(x$samples)-1,"\n")
cat("spectral statistic used:",x$statistic,"\n\n")

}
