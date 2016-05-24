###############################################################################
## Some simple tests
###############################################################################

library(RobLox)

## sample
x <- rnorm(10000, mean = -2, sd = 3)

## location and scale, radius unknown
res1 <- roblox(x, returnIC = TRUE)
checkIC(pIC(res1))

res11 <- roblox(x, returnIC = TRUE, k = 2)
checkIC(pIC(res11))

res12 <- roblox(x, returnIC = TRUE, k = 5)
checkIC(pIC(res12))

roblox(x)
roblox(x, k = 2)
roblox(x, k = 5)


## location and scale, radius interval
res2 <- roblox(x, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE)
checkIC(pIC(res2))

res21 <- roblox(x, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 2)
checkIC(pIC(res21))

res22 <- roblox(x, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 4)
checkIC(pIC(res22))

roblox(x, eps.lower = 0.15, eps.upper = 0.3)
roblox(x, eps.lower = 0.15, eps.upper = 0.3, k = 2)
roblox(x, eps.lower = 0.15, eps.upper = 0.3, k = 4)


## scale, radius interval
res3 <- roblox(x, mean = -2, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE)
checkIC(pIC(res3))

res31 <- roblox(x, mean = -2, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 3)
checkIC(pIC(res31))

res32 <- roblox(x, mean = -2, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 6)
checkIC(pIC(res32))

roblox(x, mean = -2, eps.lower = 0.15, eps.upper = 0.3)
roblox(x, mean = -2, eps.lower = 0.15, eps.upper = 0.3, k = 3)
roblox(x, mean = -2, eps.lower = 0.15, eps.upper = 0.3, k = 6)


## location, radius interval
res4 <- roblox(x, sd = 3, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE)
checkIC(pIC(res4))

res41 <- roblox(x, sd = 3, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 2)
checkIC(pIC(res41))

res42 <- roblox(x, sd = 3, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 5)
checkIC(pIC(res42))

roblox(x, sd = 3, eps.lower = 0.15, eps.upper = 0.3)
roblox(x, sd = 3, eps.lower = 0.15, eps.upper = 0.3, k = 2)
roblox(x, sd = 3, eps.lower = 0.15, eps.upper = 0.3, k = 5)


## some timings
system.time(for(i in 1:100) roblox(x, eps = 0.02))
system.time(for(i in 1:100) roblox(x))
system.time(for(i in 1:100) roblox(x, k = 2))
system.time(for(i in 1:100) roblox(x, k = 5))
system.time(for(i in 1:100) roblox(x, returnIC = TRUE))


## Samples
X <- matrix(rnorm(200, mean = -2, sd = 3), nrow = 2)
X1 <- rbind(rnorm(100, mean = -2, sd = 3), rnorm(100, mean = -1, sd = 4))

## location: rowRoblox
apply(X, 1, roblox, sd = 3)
rowRoblox(X, sd = 3)
apply(X, 1, roblox, eps = 0.06, sd = 3)
rowRoblox(X, eps = 0.06, sd = 3)

apply(X, 1, roblox, sd = 3, k = 3)
rowRoblox(X, sd = c(3, 3), k = 3)
apply(X, 1, roblox, eps = 0.06, sd = 3, k = 3)
rowRoblox(X, eps = 0.06, sd = 3, k = 3)

roblox(X1[1,], sd = 3)
roblox(X1[2,], sd = 4)
rowRoblox(X1, sd = c(3, 4))
roblox(X1[1,], eps = 0.06, sd = 3)
roblox(X1[2,], eps = 0.06, sd = 4)
rowRoblox(X1, eps = 0.06, sd = c(3, 4))

roblox(X1[1,], sd = 3, k = 4)
roblox(X1[2,], sd = 4, k = 4)
rowRoblox(X1, sd = c(3, 4), k = 4)
roblox(X1[1,], eps = 0.06, sd = 3, k = 4)
roblox(X1[2,], eps = 0.06, sd = 4, k = 4)
rowRoblox(X1, eps = 0.06, sd = c(3, 4), k = 4)


## scale: rowRoblox
apply(X, 1, roblox, mean = -2)
rowRoblox(X, mean = -2)
apply(X, 1, roblox, eps = 0.049, mean = -2)
rowRoblox(X, eps = 0.049, mean = -2)

apply(X, 1, roblox, mean = -2, k = 3)
rowRoblox(X, mean = c(-2, -2), k = 3)
apply(X, 1, roblox, eps = 0.049, mean = -2, k = 3)
rowRoblox(X, eps = 0.049, mean = c(-2, -2), k = 3)

roblox(X1[1,], mean = -2)
roblox(X1[2,], mean = -1)
rowRoblox(X1, mean = c(-2, -1))
roblox(X1[1,], eps = 0.049, mean = -2)
roblox(X1[2,], eps = 0.049, mean = -1)
rowRoblox(X1, eps = 0.049, mean = c(-2, -1))

roblox(X1[1,], mean = -2, k = 4)
roblox(X1[2,], mean = -1, k = 4)
rowRoblox(X1, mean = c(-2, -1), k = 4)
roblox(X1[1,], eps = 0.049, mean = -2, k = 4)
roblox(X1[2,], eps = 0.049, mean = -1, k = 4)
rowRoblox(X1, eps = 0.049, mean = c(-2, -1), k = 4)

## location and scale: rowRoblox
apply(X, 1, roblox)
rowRoblox(X)
apply(X, 1, roblox, eps = 0.057)
rowRoblox(X, eps = 0.057)

apply(X, 1, roblox, k = 3)
rowRoblox(X, k = 3)
apply(X, 1, roblox, eps = 0.057, k = 3)
rowRoblox(X, eps = 0.057, k = 3)


## Samples
X <- t(X)
X1 <- t(X1)

## location: colRoblox
apply(X, 2, roblox, sd = 3)
colRoblox(X, sd = 3)
apply(X, 2, roblox, eps = 0.06, sd = 3)
colRoblox(X, eps = 0.06, sd = 3)

apply(X, 2, roblox, sd = 3, k = 3)
colRoblox(X, sd = c(3, 3), k = 3)
apply(X, 2, roblox, eps = 0.06, sd = 3, k = 3)
colRoblox(X, eps = 0.06, sd = 3, k = 3)

roblox(X1[,1], sd = 3)
roblox(X1[,2], sd = 4)
colRoblox(X1, sd = c(3, 4))
roblox(X1[,1], eps = 0.06, sd = 3)
roblox(X1[,2], eps = 0.06, sd = 4)
colRoblox(X1, eps = 0.06, sd = c(3, 4))

roblox(X1[,1], sd = 3, k = 4)
roblox(X1[,2], sd = 4, k = 4)
colRoblox(X1, sd = c(3, 4), k = 4)
roblox(X1[,1], eps = 0.06, sd = 3, k = 4)
roblox(X1[,2], eps = 0.06, sd = 4, k = 4)
colRoblox(X1, eps = 0.06, sd = c(3, 4), k = 4)


## scale: colRoblox
apply(X, 2, roblox, mean = -2)
colRoblox(X, mean = -2)
apply(X, 2, roblox, eps = 0.049, mean = -2)
colRoblox(X, eps = 0.049, mean = -2)

apply(X, 2, roblox, mean = -2, k = 3)
colRoblox(X, mean = c(-2, -2), k = 3)
apply(X, 2, roblox, eps = 0.049, mean = -2, k = 3)
colRoblox(X, eps = 0.049, mean = c(-2, -2), k = 3)

roblox(X1[,1], mean = -2)
roblox(X1[,2], mean = -1)
colRoblox(X1, mean = c(-2, -1))
roblox(X1[,1], eps = 0.049, mean = -2)
roblox(X1[,2], eps = 0.049, mean = -1)
colRoblox(X1, eps = 0.049, mean = c(-2, -1))

roblox(X1[,1], mean = -2, k = 4)
roblox(X1[,2], mean = -1, k = 4)
colRoblox(X1, mean = c(-2, -1), k = 4)
roblox(X1[,1], eps = 0.049, mean = -2, k = 4)
roblox(X1[,2], eps = 0.049, mean = -1, k = 4)
colRoblox(X1, eps = 0.049, mean = c(-2, -1), k = 4)

## location and scale: colRoblox
apply(X, 2, roblox)
colRoblox(X)
apply(X, 2, roblox, eps = 0.057)
colRoblox(X, eps = 0.057)

apply(X, 2, roblox, k = 3)
colRoblox(X, k = 3)
apply(X, 2, roblox, eps = 0.057, k = 3)
colRoblox(X, eps = 0.057, k = 3)


## some timings
X <- matrix(rnorm(1e5, mean = -1, sd = 3), ncol = 100)
system.time(apply(X, 1, roblox, eps = 0.02))
## uses rowMedians of package Biobase if available
system.time(rowRoblox(X, eps = 0.02))

system.time(apply(X, 1, roblox))
## uses rowMedians of package Biobase if available
system.time(rowRoblox(X))

M <- apply(X, 1, median)
S <- apply(X, 1, mad)
init <- cbind(M, S)
system.time(rowRoblox(X, initial.est = init))
