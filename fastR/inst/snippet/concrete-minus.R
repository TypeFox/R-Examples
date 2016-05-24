# modify data by dropping first observation
concretemod <- concrete[-1,]
concrete.lmmod <- lm(strength ~ limestone + water, concretemod)
coef(concrete.lmmod)
y <- concretemod$strength
n <- length(y); v0 <- rep(1,n)
v1 <- concretemod$limestone - mean(concretemod$limestone)
v2 <- concretemod$water - mean(concretemod$water)
project(y,v0,type='v')
mean(y)
project(y,v1,type='l') / vlength(v1)
project(y,v2,type='l') / vlength(v2)
ef0 <- project(y,v0,type='v')
ef1 <- project(y,v1,type='v')
ef2 <- project(y,v2,type='v')
ef0 + ef1 + ef2
predict(concrete.lmmod)
