RSVc <-
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
m_t<-mean(t)
tord<-sort(t)
for (a in 2:H)
{
tord[a]<-tord[a-1]+tord[a]
}
tord<-tord/(m_t*H)
xaxis<-cumsum(p)/N
xaxis
RSV<-c(0,tord)
p<-c(0,xaxis)
RSVc<-list(p,RSV)
names(RSVc)<-c("RSV (maximum inequality) x-axis points","RSV (maximum inequality) y-axis points")
x<-RSVc$'RSV (maximum inequality) x-axis points'
y<-RSVc$'RSV (maximum inequality) y-axis points'
if (plot)
plot(x,y,col="black",xlab="p",ylab="L(p)",main="RSV curve of maximum inequality",xaxs="i",type="l")
Last<-length(RSVc$'RSV (maximum inequality) x-axis points')
segments(0,0,1,0,col="black")
segments(RSVc$'RSV (maximum inequality) x-axis points'[1],RSVc$'RSV (maximum inequality) y-axis points'[1],RSVc$'RSV (maximum inequality) x-axis points'[Last],RSVc$'RSV (maximum inequality) x-axis points'[Last],col="black",lwd="1")
for (a in 1:(Last-1))
{
segments(RSVc$'RSV (maximum inequality) x-axis points'[a],RSVc$'RSV (maximum inequality) y-axis points'[a],RSVc$'RSV (maximum inequality) x-axis points'[a+1],RSVc$'RSV (maximum inequality) y-axis points'[a+1],col="black",lwd="2")
}
RSVc
}
