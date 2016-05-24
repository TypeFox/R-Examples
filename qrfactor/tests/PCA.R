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

