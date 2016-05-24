y <- concretemod$strength
# make fits using v1 and w2
ef0 <- project(y,v0,type='v')
ef1 <- project(y,v1,type='v')
ef2 <- project(y,w2,type='v')
ef0 + ef1 + ef2
# now try w1 and v2
ef0 <- project(y,v0,type='v')
ef1 <- project(y,w1,type='v')
ef2 <- project(y,v2,type='v')
ef0 + ef1 + ef2
# should match what lm() produces
predict(concrete.lmmod)
