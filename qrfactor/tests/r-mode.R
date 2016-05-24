library(qrfactor)
#data must be numeric, avoid using categorical and characters
data(UScereal, package="MASS")
variables=c("calories","protein","sodium","carbo","sugars","potassium")
data=UScereal[variables]

#R-mode Factor Analysis
cat("qrfactor R mode Loadings\n")
mod=qrfactor(data)
mod$r.loading
plot(mod,type="loadings",plot="r")
plot(mod,type="loadings",plot="all")

