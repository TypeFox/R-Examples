best.unbal.haar.bu <-
function(x, stretch = length(x)) {

n <- length(x)
weights <- rep(1, n)
detail <- matrix(0, 3, n-1)

for (i in 1:(n-1)) {

tmp <- x[i:min(i+stretch-1,n)]

weights.tmp <- weights[i:min(i+stretch-1,n)]

m <- length(tmp)

a <- weights.tmp[1:(m-1)]
b <- weights.tmp[2:m]
h1 <- (1 + a^2 / b^2)^(-1/2)
h2 <- -(1 + b^2 / a^2)^(-1/2)

l1 <- -h2
l2 <- h1

tmp.dif <- h1 * tmp[1:(m-1)] + h2 * tmp[2:m]
tmp.dif.min <- min(which(abs(tmp.dif) == min(abs(tmp.dif))))

detail[1,i] <- tmp.dif.min
detail[2,i] <- l1[tmp.dif.min]
detail[3,i] <- tmp.dif[tmp.dif.min]

tmp.mod <- tmp[-tmp.dif.min]

tmp.mod[tmp.dif.min] <- l1[tmp.dif.min] * tmp[tmp.dif.min] + l2[tmp.dif.min] * tmp[tmp.dif.min+1]
x[(i+1):min(i+stretch-1,n)] <- tmp.mod

weights.tmp.mod <- weights.tmp[-tmp.dif.min]
weights.tmp.mod[tmp.dif.min] <- l1[tmp.dif.min] * weights.tmp[tmp.dif.min] + l2[tmp.dif.min] * weights.tmp[tmp.dif.min+1]
weights[(i+1):min(i+stretch-1,n)] <- weights.tmp.mod

}

smooth <- x[n]

return(list(smooth = smooth, detail = detail))

}
