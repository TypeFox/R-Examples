library(qrfactor)
#data must be numeric, avoid using categorical and characters
data(UScereal, package="MASS")
variables=c("calories","protein","sodium","carbo","sugars","potassium")
data=UScereal[variables]

#Q-mode Factor Analysis
mod=qrfactor(data)
mod$q.loading
summary(mod)
plot(mod,type="loadings",plot="q")
plot(mod,type="scores",plot="q")

