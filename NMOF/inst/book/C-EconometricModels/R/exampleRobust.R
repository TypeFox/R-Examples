# exampleRobust.R -- version 2010-12-29
# ... continues 'comparisonLMS.R'

n <- 100L  # number of observations
p <- 10L   # number of regressors
constant <- TRUE; sigma <- 5; oFrac  <- 0.15
h <- 70L   # ... or use something like floor((n+1)/2)
aux <- createData(n,p,constant,sigma,oFrac) 
X <- aux$X; y <- aux$y

trials <- 100L
res1 <- numeric(trials)
for (t in 1L:trials){
    modl <- lqs(y ~ X[ ,-1L], adjust = TRUE, 
    nsamp = "best", method = "lqs", quantile = h)
    res1[t] <- sort((y - X %*% as.matrix(coef(modl)))^2)[h]
}
#
res2 <- numeric(trials)
for (t in 1L:trials){
    modl <- lqs(y ~ X[ ,-1L], adjust = TRUE, 
    nsamp = 10000L, method = 'lqs', quantile = h)
    res2[t] <- sort((y - X %*% as.matrix(coef(modl)))^2)[h]
}
#
res3 <- numeric(trials)
for (t in 1L:trials){
    modl <- lqs(y ~ X[ ,-1L], adjust = TRUE, 
    nsamp = 100000L, method = 'lqs', quantile = h)
    res3[t] <- sort((y - X %*% as.matrix(coef(modl)))^2)[h]
}
xx <- pretty(res1, res2, res3)
Y <- sort(res1); N <- trials
plot( c(Y[1L], Y), (0L:N)/N, type = 's', col = grey(0.6), 
    xlim = c(min(xx),max(xx)))
Y <- sort(res2); N <- trials
lines(c(Y[1L] ,Y), (0L:N)/N, type = 's', col = grey(0.2))
Y <- sort(res3); N <- trials
lines(c(Y[1L] ,Y), (0L:N)/N, type = 's', col = grey(0.0))
