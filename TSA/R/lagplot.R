`lagplot` <-
function(x,lag.max=6,deg=1,nn=0.7, method=c("locfit","gam","both")[1]){
require(locfit)
require(mgcv)
par(mfrow=c(3,2),mar=c(3,2,3,2))
if(method=='gam'){
for (k in 1:lag.max){
g1=gam(x~zlag(x,k))
pg1=predict(g1,type='response')
plot(y=x,x=zlag(x,k),xlab='',ylab='',asp=1,main=paste(paste('lag',k,sep='-'),
"regression plot"))
lines(y=pg1,x=na.omit(zlag(x,k)),lty=2,col='red')
} } else {
for (k in 1:lag.max){
l1=locfit(x~lp(zlag(x,k),deg=deg,nn=nn))
plot(l1,get.data=TRUE,asp=1,main=paste(paste('lag',k,sep='-'),"regression plot"))
if(method=="both") {
g1=gam(x~zlag(x,k))
pg1=predict(g1,type='response')
lines(y=pg1,x=na.omit(zlag(x,k)),lty=2,col='red')

}
}
}
invisible()
}

