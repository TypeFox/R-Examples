if (interactive()) par.default <- par(ask=TRUE)

#-----  Test Gaussian Distribution -----
#             fixed variance

Pi <- matrix(c(1/2, 1/2,   0,
               1/3, 1/3, 1/3,
                 0, 1/2, 1/2),
             byrow=TRUE, nrow=3)
p1 <- c(1, 6, 3)
p2 <- c(0.5, 1, 0.5)
n <- 1000

pm <- list(mean=p1)
pn <- list(sd=rep(0.5, n))

n <- 1000
x <- dthmm(NULL, Pi, c(0,1,0), "norm", pm=pm, pn=pn)

x <- simulate(x, nsim=n)

#    use above parameter values as initial values
y <- BaumWelch(x)

#   check parameter estimates
print(summary(y))
print(sum(y$delta))
print(y$Pi %*% rep(1, ncol(y$Pi)))


#-----  Test Gaussian Distribution -----
#   residual process, zero mean, Markov driven variance

Pi <- matrix(c(0.9, 0.1,
               0.9, 0.1),
             byrow=TRUE, nrow=2)
n <- 1000

pn <- list(mean=rep(0, n))
pm <- list(sd=c(1, 3))

n <- 1000
x <- dthmm(NULL, Pi, c(0,1), "norm", pm=pm, pn=pn)

x <- simulate(x, nsim=n)

#   check transition probabilities
m <- nrow(Pi)
estPi <- table(x$y[-n], x$y[-1])
rowtotal <- estPi %*% matrix(1, nrow=m, ncol=1)
estPi <- diag(as.vector(1/rowtotal)) %*% estPi
print(estPi)

#    use above parameter values as initial values
y <- BaumWelch(x)

#   check parameter estimates
print(summary(y))
print(sum(y$delta))
print(y$Pi %*% rep(1, ncol(y$Pi)))

plot(1:n, x$x, type="l", xlab="Time", ylab="")
points((1:n)[x$y==2], x$x[x$y==2], pch=15, col="red")
title(sub="Those marked in red are in state 2")

states <- Viterbi(y)
wrong <- (states!=x$y)
print(table(x$y, wrong))

par(mfrow=c(1, 2))
s1 <- hist(x$x[x$y==1], col="gray70", xlab="",
           main="Observed X Values When in State 1")
hist(x$x[x$y==1 & wrong], breaks=s1$breaks, add=TRUE, col="red")
title(sub="Red denotes the number of incorrectly determined states by Viterbi")
box()
s2 <- hist(x$x[x$y==2], col="gray70", xlab="",
           main="Observed X Values When in State 2")
hist(x$x[x$y==2 & wrong], breaks=s2$breaks, add=TRUE, col="red")
title(sub="Red denotes the number of incorrectly determined states by Viterbi")
box()
par(mfrow=c(1, 1))

if (interactive()) par(par.default)
