library(qrfactor)
#data must be numeric, avoid using categorical and characters
data(UScereal, package="MASS")
variables=c("calories","protein","sodium","carbo","sugars","potassium")
data=UScereal[variables]
#create object with observation number starting with nothing. 
mod1 <- qrfactor(data)
#print object
mod1
#extract the loadings
loadings=mod1$loadings
loadings
mod1$r.loading
mod1$q.loading

#summary of the loadings
summary(mod1)
#plot the first two axes of the loadings
plot(mod1)

#plotting other axes for mod1 for axes 3 and 4
plot(mod1,factors=c(3,4))
plot(mod1,factors=c(1,3))
#plotting scores fro axes 2 and 3
plot(mod1,factors=c(3,4),type='scores')

#Having multiple plots on page
plot(mod1,factors=c(3,4),plot="all")
plot(mod1,factors=c(1,3),plot="qr")
#plotting scores fro axes 2 and 3
plot(mod1,factors=c(3,4),type='scores',plot="r")
plot(mod1)

#principal coordinate analysis:
plot(qrfactor(data),plot="all",type='coord',factors=c(1,3))
plot(qrfactor(data),type='coord',factors=c(2,3))


#Multi dimensional scaling can be simulated as: 
plot(qrfactor(data,scale='n'),plot="r",type="mds",factors=c(1,2))
plot(qrfactor(data,scale='n'),plot="all",type="mds",factors=c(3,2))



