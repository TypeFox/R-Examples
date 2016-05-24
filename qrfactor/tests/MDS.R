library(qrfactor)
#data must be numeric, avoid using categorical and characters
data(UScereal, package="MASS")
variables=c("calories","protein","sodium","carbo","sugars","potassium")
data=UScereal[variables]
#Multi dimensional scaling can be simulated as: 
plot(qrfactor(data,scale='n'),plot="r",type="mds",factors=c(1,2))
plot(qrfactor(data,scale='n'),plot="all",type="mds",factors=c(3,2))

