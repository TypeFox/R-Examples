GRSVc <-
function(z,p=rep(1,length(z)),plot=FALSE){
zp<-z*p
N<-sum(p)
H<-length(z)
t_neg<-subset(zp,zp<=0)
T_neg<-sum(t_neg)
T_neg
T_pos<-sum(zp)+abs(T_neg)
T_pos
t_intra<-rep(0,H-2)
t<-c(T_neg,t_intra,T_pos)
tord<-sort(t)
for (a in 2:H)
{
tord[a]<-tord[a-1]+tord[a]
}
tord<-tord/N
xaxis<-cumsum(p)/N
xaxis
RSV<-c(0,tord)
p<-c(0,xaxis)
GRSVc<-list(p,RSV)
names(GRSVc)<-c("Generalized RSV (maximum inequality) x-axis points","Generalized RSV (maximum inequality) y-axis points")
x<-GRSVc$'Generalized RSV (maximum inequality) x-axis points'
y<-GRSVc$'Generalized RSV (maximum inequality) y-axis points'
if (plot)
plot(x,y,col="black",xlab="p",ylab="L(p)",lwd="2",main="Generalized RSV curve of maximum inequality",xaxs="i",type="l")
Last<-length(GRSVc$'Generalized RSV (maximum inequality) x-axis points')
segments(0,0,1,0,col="black")
segments(GRSVc$'Generalized RSV (maximum inequality) x-axis points'[1],GRSVc$'Generalized RSV (maximum inequality) y-axis points'[1],GRSVc$'Generalized RSV (maximum inequality) x-axis points'[Last],GRSVc$'Generalized RSV (maximum inequality) y-axis points'[Last],col="black",lwd="1")
for (a in 1:(Last-1))
{
segments(GRSVc$'Generalized RSV (maximum inequality) x-axis points'[a],GRSVc$'Generalized RSV (maximum inequality) y-axis points'[a],GRSVc$'Generalized RSV (maximum inequality) x-axis points'[a+1],GRSVc$'Generalized RSV (maximum inequality) y-axis points'[a+1],col="black",lwd="2")
}
GRSVc
}
