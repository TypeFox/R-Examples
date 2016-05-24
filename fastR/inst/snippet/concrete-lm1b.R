y <- concrete$strength
ef0 <- project(y,v0,type='v')
ef1 <- project(y,v1,type='v')
ef2 <- project(y,v2,type='v')
ef0 + ef1 + ef2
predict(concrete.lm1)
fitted(concrete.lm1)
