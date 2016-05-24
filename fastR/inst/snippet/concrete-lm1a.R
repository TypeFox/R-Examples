y <- concrete$strength
n <- length(y); v0 <- rep(1,n)
v1 <- concrete$limestone - mean(concrete$limestone)
v2 <- concrete$water - mean(concrete$water)
project(y,v0,type='v')
mean(y)
project(y,v1,type='l') / vlength(v1)
project(y,v2,type='l') / vlength(v2)
