hindcast<-function(Dhigh,hightimes,Dlow,lowtimes,predtimes,filter=c(1),ebeta,vbeta,SARIMA){

#hightimes<-1:16
#lowtimes<-16+1:5*2
#predtimes<-16+1:11

# Gather series lengths
deltat<-abs(lowtimes[2]-lowtimes[1])
nhigh<-length(hightimes)
nlow<-length(lowtimes)
nfilter<-length(filter)
npred<-length(predtimes)
# Identify hidden past
lowtimes.all<-(min(lowtimes)-nfilter+1):max(lowtimes)
nlow.all<-length(lowtimes.all)

# Compute autocovariances
timerange<-range(lowtimes.all,hightimes,predtimes)
maxlag<-timerange[2]-timerange[1]
k<-logspec2cov(ebeta,vbeta,SARIMA=SARIMA,lags=maxlag)
D<-as.matrix(dist(c(lowtimes.all,hightimes,predtimes)))
K<-matrix(k[D+1],nrow(D),nrow(D))

# Construct observation matrix
B<-t(matrix(c(rep(c(rev(filter),rep(0,nlow.all-nfilter+deltat)),nlow-1),rev(filter)),nlow.all,nlow))

# Construct variance matrix
V<-rbind(cbind(B%*%K[1:nlow.all,1:nlow.all]%*%t(B),B%*%K[1:nlow.all,nlow.all+1:(nhigh+npred)]),
cbind(t(B%*%K[1:nlow.all,nlow.all+1:(nhigh+npred)]),K[nlow.all+1:(nhigh+npred),nlow.all+1:(nhigh+npred)]))

# Take bits out
CBD<-V[nlow+nhigh+1:npred,1:(nlow+nhigh)]
VD<-V[1:(nlow+nhigh),1:(nlow+nhigh)]

# Compute hindcast
VDinv<-solve(VD)
hindcast<-CBD%*%VDinv%*%c(Dlow,Dhigh)
var.hindcast<-V[nlow+nhigh+1:npred,nlow+nhigh+1:npred]-CBD%*%VDinv%*%t(CBD)

list(hindcast=hindcast,var.hindcast=var.hindcast)
}
