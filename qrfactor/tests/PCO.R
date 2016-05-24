library(qrfactor)
#data must be numeric, avoid using categorical and characters
data(UScereal, package="MASS")
variables=c("calories","protein","sodium","carbo","sugars","potassium")
data=UScereal[variables]

#principal coordinate analysis:
plot(qrfactor(data),plot="all",type='coord',factors=c(1,3))
plot(qrfactor(data),type='coord',factors=c(2,3))

