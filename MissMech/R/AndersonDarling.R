AndersonDarling <- function(data, number.cases) { 
# Not adjusted for ties
# x is data vector
# ni is number of cases in each group
  x <- data
  ni <- number.cases
 if(length(ni)<2)
 {
     cat("Warning: Not enough groups for AndersonDarling test.")
     stop("")
 }
 
 k <- length(ni)
 ni.z <- c(0, cumsum(ni))
 n <- length(x)
 x.sort <- sort(x)
 x.sort <- x.sort[1:(n-1)]
 ind <- which(duplicated(x.sort) == 0)
 counts <- c(ind, length(x.sort) + 1) - c(0, ind)
 hj <- counts[2 : (length(ind) + 1)]
 hn <- cumsum(hj)
 zj <- x.sort[ind]
 adk <- 0
 adk.all <- matrix(0, k, 1) # to keep contribution of kth group
 for(i in 1:k)
 {
    ind <- (ni.z[i] + 1) : ni.z[i + 1]
    templist <- expand.grid(zj, x[ind])
    b <- templist[, 1] == templist[, 2]
    fij <- apply(matrix(b, length(zj)), 1, sum)
    mij <- cumsum(fij)
    num <- (n * mij - ni[i] * hn) ^ 2
    dem <- hn*(n - hn)
    adk.all[i] <- (1 / ni[i] * sum(hj * (num / dem)))
    adk <- adk + adk.all[i]
 }
 adk <- (1 / n) * adk
 adk.all <- adk.all / n
#Exact sample variance of the k-sample Anderson-Darling
#Finding Variance of the statistics
 j <- sum(1 / ni)
 i <- seq(1:(n - 1))
 h <- sum(1 / i)
 g <- 0
 for (i in 1:(n - 2)) {
     g <- g + (1 / (n - i)) * sum(1 / seq((i + 1), (n - 1)))
 }
 a <- (4 * g - 6) * (k - 1) + (10 - 6 * g) * j
 b <- (2 * g - 4) * k ^ 2 + 8 * h * k + 
      (2 * g - 14 * h - 4) * j - 8 * h + 4 * g - 6
 c <- (6 * h + 2 * g - 2) * k ^ 2 + (4 * h - 4 * g + 6) * k +
      (2 * h - 6) * j + 4 * h
 d <- (2 * h + 6) * k ^ 2 - 4 * h * k
 var.adk <- ((a * n ^ 3) + (b * n ^ 2) + (c * n) + d) /
            ((n - 1) * (n - 2) * (n - 3))
 if(var.adk<0) var.adk=0
 adk.s <- (adk - (k - 1)) / sqrt(var.adk) 
#k-sample Anderson-Darling P-value calculation by an extrapolate-interpolate
#procedure
 a0 <- c(0.25, 0.10, 0.05, 0.025, 0.01)
 b0 <- c(0.675, 1.281, 1.645, 1.96, 2.326)
 b1 <- c(-0.245, 0.25, 0.678 ,1.149, 1.822)
 b2 <- c(-0.105, -0.305, -0.362, -0.391, -0.396)
 c0 <- log((1 - a0) / a0)
 qnt <- b0 + b1 / sqrt(k - 1) + b2 / (k - 1)
 if (adk.s <= qnt[3]) {
    ind <- seq(1:4)
 } else {
    ind <- seq(2:5)
 }
 yy <- spline(qnt[ind], c0[ind], xout = adk.s)$y
 p <- 1 / (1 + exp(yy))
 list(pn = p, adk.all = adk.all, adk = adk, var.sdk = var.adk)
}#end function

