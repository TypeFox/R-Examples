POTevents.fun <-
function(T,thres, date=NULL)
{
if (is.null(date)) date<-c(1:length(T))
date<-as.matrix(date)
if ((is.null(date)==FALSE)&(dim(date)[1]!=length(T)) ) 
stop('T and date must have the same number of observations')
exc<-(T>thres)
inrachtx<-(c(0,diff(exc))==1)

Pi<-c(1:length(T))[inrachtx==1]
numerachtx<-cumsum(inrachtx)[exc==1]

intentx<-(T-thres)[exc==1]

Im<-tapply(intentx,INDEX=numerachtx, FUN=mean)
Ix<-tapply(intentx,INDEX=numerachtx, FUN=max)
L<-tapply(intentx,INDEX=numerachtx, FUN=length)
Px<-Pi+tapply(intentx,INDEX=numerachtx, FUN=which.max)-1

inddat<-1-exc
inddat[Px]<-1

datePi<-date[Pi,]
datePx<-date[Px,]

cat('Number of events: ',length(Im),fill=TRUE)
cat('Number of excesses over threshold', thres,':',sum(exc),fill=TRUE)

return(list(Pi=Pi, datePi=datePi, Px=Px, datePx=datePx, Im=Im,Ix=Ix,L=L,inddat=inddat,T=T,thres=thres, date=date))

}
