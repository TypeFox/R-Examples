# using v1 gives coefficient in model with only limestone as a predictor
project(y,v1,type='l') / vlength(v1)
coef(lm(strength ~ limestone, concretemod))
# using v2 gives coefficient in model with only water as a predictor
project(y,v2,type='l') / vlength(v2)
coef(lm(strength ~ water, concretemod))
# using w1 and w2 gives coefficients in the additive model 
project(y,w1,type='l') / vlength(w1)
project(y,w2,type='l') / vlength(w2)
coef(concrete.lmmod)
