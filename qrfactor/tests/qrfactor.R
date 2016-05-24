library(qrfactor)
#data must be numeric, avoid using categorical and characters
data(UScereal, package="MASS")
variables=c("calories","protein","sodium","carbo","sugars","potassium")
data=UScereal[variables]


#PCA
mod=qrfactor(data)
cat("PCA standard deviation for qrfactor\n")
sqrt(mod$eigen.value)
cat("PCA Loadings for qrfactor\n")
mod$pca
#plot PCA
plot(mod,type="PCA",plot="all")

#R-mode Factor Analysis
cat("qrfactor R mode Loadings\n")
mod=qrfactor(data)
mod$r.loading
plot(mod,type="loadings",plot="r")
plot(mod,type="loadings",plot="all")

#Q-mode Factor Analysis
mod=qrfactor(data)
mod$q.loading
summary(mod)
plot(mod,type="loadings",plot="q")
plot(mod,type="scores",plot="q")

#Simulataneous R and Q-mode Factor Analysis
mod=qrfactor(data)
mod$loadings
plot(mod,type="loadings")
plot(mod,type="scores")
plot(mod,type="loadings",plot="all")
plot(mod,type="scores",plot="qr")

#principal coordinate analysis:
plot(qrfactor(data),plot="all",type='coord',factors=c(1,3))
plot(qrfactor(data),type='coord',factors=c(2,3))


#Multi dimensional scaling can be simulated as: 
plot(qrfactor(data,scale='n'),plot="r",type="mds",factors=c(1,2))
plot(qrfactor(data,scale='n'),plot="all",type="mds",factors=c(3,2))

