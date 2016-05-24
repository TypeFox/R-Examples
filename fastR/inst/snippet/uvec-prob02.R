x <- c(3,4,5,8)
mean(x)
var(x)
l <- c(); P<- c()
for (i in 1:4) {
    l[i] <- as.numeric(x %*% uvec(i,4))      # proj coefficient
    P[[i]] <- l[i] * uvec(i,4)               # proj vector 
}
l                                            # proj coefficients
P                                            # proj vectors
# next two should be the same value
l[2]^2 + l[3]^2 + l[4]^2                          
3 * var(x)
