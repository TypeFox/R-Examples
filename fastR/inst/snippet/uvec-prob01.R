x <- c(3,4,4,7,7)
mean(x)
v <- x - mean(x)
lx <- c(); lv <- c()
for (i in 1:5) {
    lx[i] <- as.numeric(x %*% uvec(i,5))     # proj coefficient
    lv[i] <- as.numeric(v %*% uvec(i,5))     # proj coefficient
}
# all but first should match exactly; lv[1] == 0
lx                                           # projecting from x
lv                                           # projecting from v
