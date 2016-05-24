x <- c(1:20,20)
boxplot.stats(x)
iqr <- 16 - 6
# not an outlier
x <- c(1:20, 16 + 1.5 * iqr)
boxplot.stats(x)$out
# now it is an outlier
x <- c(1:20, 16 + 1.51 * iqr)
boxplot.stats(x)$out
###
