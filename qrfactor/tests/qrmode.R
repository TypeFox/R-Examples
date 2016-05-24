library(qrfactor)
#data must be numeric, avoid using categorical and characters
data(UScereal, package="MASS")
variables=c("calories","protein","sodium","carbo","sugars","potassium")
data=UScereal[variables]

#Simulataneous R and Q-mode Factor Analysis
mod=qrfactor(data)
mod$loadings
plot(mod,type="loadings")
plot(mod,type="scores")
plot(mod,type="loadings",plot="all")
plot(mod,type="scores",plot="qr")

