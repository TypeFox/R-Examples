pl <-
function(BOR1,BTR1,BOR2,BTR2,MainLabel1="",MainLabel2="",HA1="",
HA2="",VA1="",VA2="")
{
par(mfrow=c(1,2))
plot(BOR1,BTR1,type="n",main=MainLabel1,xlab=HA1,ylab=VA1,font=15)
abline(h=mean(BTR1, na.rm=TRUE),col=2,lwd=2)
abline(v=mean(BOR1, na.rm=TRUE),col=4,lwd=3)

(a=data.frame(BOR1,BTR1))
(m1=a[a[,1]<mean(BOR1),]);(m2=m1[m1[,2]<mean(BTR1),])
points(m2,pch=1,lwd=2,cex=3)


(m1=a[a[,1]<mean(BOR1),]);(m2=m1[m1[,2]>mean(BTR1),])
points(m2,pch=1,lwd=2,cex=3)


(m1=a[a[,1]>mean(BOR1),]);(m2=m1[m1[,2]>mean(BTR1),])
points(m2,pch=1,lwd=2,cex=3)


(m1=a[a[,1]>mean(BOR1),]);(m2=m1[m1[,2]<mean(BTR1),])
points(m2,pch=1,lwd=2,cex=3)
text(BOR1,BTR1,1:length(BOR1))

plot(BOR2,BTR2,type="n",main=MainLabel2,xlab=HA2,ylab=VA2,font=15)
abline(h=mean(BTR2, na.rm=TRUE),col=2,lwd=2)
abline(v=mean(BOR2, na.rm=TRUE),col=4,lwd=3)
abline(h=mean(BTR1, na.rm=TRUE),col=2,lty=2,lwd=2)
abline(v=mean(BOR1, na.rm=TRUE),lty=3,col=4,lwd=3)
#Remainding in Zone I
(a=data.frame(BOR2,BTR2,BOR1,BTR1))
n1=a[a[,1]<mean(BOR2, na.rm=TRUE)&a[,3]<mean(BOR1, na.rm=TRUE),]
n2=n1[n1[ ,2]<mean(BTR2, na.rm=TRUE)&n1[ ,4]<mean(BTR1, na.rm=TRUE),]
points(n2,pch=1,lwd=2,cex=3)

#Transition from Zone II to Zone I
(n1=a[a[,1]<mean(BOR2, na.rm=TRUE)&a[,3]<mean(BOR1, na.rm=TRUE),])
(n2=n1[n1[ ,2]<mean(BTR2, na.rm=TRUE)&n1[ ,4]>mean(BTR1, na.rm=TRUE),])
points(n2,pch=19,col="red",lwd=2,cex=3)

#Transition from Zone III to Zone I
(n1=a[a[,1]<mean(BOR2, na.rm=TRUE)&a[,3]>mean(BOR1, na.rm=TRUE),])
(n2=n1[n1[ ,2]<mean(BTR2, na.rm=TRUE)&n1[ ,4]>mean(BTR1, na.rm=TRUE),])
points(n2,pch=19, col="red",lwd=2,cex=3)

#Transition from Zone IV to Zone I
(n1=a[a[,1]<mean(BOR2, na.rm=TRUE)&a[,3]>mean(BOR1, na.rm=TRUE),])
(n2=n1[n1[ ,2]<mean(BTR2, na.rm=TRUE)&n1[ ,4]<mean(BTR1, na.rm=TRUE),])
points(n2,pch=19, col="red",lwd=2,cex=3)

#Transition from Zone I to Zone II
(n1=a[a[,1]<mean(BOR2, na.rm=TRUE)&a[,3]<mean(BOR1, na.rm=TRUE),])
(n2=n1[n1[ ,2]>mean(BTR2, na.rm=TRUE)&n1[ ,4]<mean(BTR1, na.rm=TRUE),])
points(n2,pch=19,col=6,lwd=2,cex=3)

#Remainding in Zone II
(n1=a[a[,1]<mean(BOR2, na.rm=TRUE)&a[,3]<mean(BOR1, na.rm=TRUE),])
(n2=n1[n1[ ,2]>mean(BTR2, na.rm=TRUE)&n1[ ,4]>mean(BTR1, na.rm=TRUE),])
points(n2,pch=1,lwd=2,cex=3)
 
#Transition from Zone III to Zone II
(n1=a[a[,1]<mean(BOR2, na.rm=TRUE)&a[,3]>mean(BOR1, na.rm=TRUE),])
(n2=n1[n1[ ,2]>mean(BTR2, na.rm=TRUE)&n1[ ,4]>mean(BTR1, na.rm=TRUE),])
points(n2,pch=19,col="greenyellow",lwd=2,cex=3)

#Transition from Zone IV to Zone II
(n1=a[a[,1]<mean(BOR2, na.rm=TRUE)&a[,3]>mean(BOR1, na.rm=TRUE),])
(n2=n1[n1[ ,2]>mean(BTR2, na.rm=TRUE)&n1[ ,4]<mean(BTR1, na.rm=TRUE),])
points(n2,pch=1,lwd=2,cex=3)

#Transition from Zone I to Zone III
(n1=a[a[,1]>mean(BOR2, na.rm=TRUE)&a[,3]<mean(BOR1, na.rm=TRUE),])
(n2=n1[n1[ ,2]>mean(BTR2, na.rm=TRUE)&n1[ ,4]<mean(BTR1, na.rm=TRUE),])
points(n2,pch=19,col=11,lwd=2,cex=3)

#Transition from Zone II to Zone III
(n1=a[a[,1]>mean(BOR2, na.rm=TRUE)&a[,3]<mean(BOR1, na.rm=TRUE),])
(n2=n1[n1[ ,2]>mean(BTR2, na.rm=TRUE)&n1[ ,4]>mean(BTR1, na.rm=TRUE),])
points(n2,pch=19,col=11,lwd=2,cex=3)

#Remainding in Zone III
(n1=a[a[,1]>mean(BOR2, na.rm=TRUE)&a[,3]>mean(BOR1, na.rm=TRUE),])
(n2=n1[n1[ ,2]>mean(BTR2, na.rm=TRUE)&n1[ ,4]>mean(BTR1, na.rm=TRUE),])
points(n2,pch=1,lwd=2,cex=3)

#Transition from Zone IV to Zone III
(n1=a[a[,1]>mean(BOR2, na.rm=TRUE)&a[,3]>mean(BOR1, na.rm=TRUE),])
(n2=n1[n1[ ,2]>mean(BTR2, na.rm=TRUE)&n1[ ,4]<mean(BTR1, na.rm=TRUE),])
points(n2,pch=19,col=11,lwd=2,cex=3)

#Transition from Zone I to Zone IV
(n1=a[a[,1]>mean(BOR2, na.rm=TRUE)&a[,3]<mean(BOR1, na.rm=TRUE),])
(n2=n1[n1[ ,2]<mean(BTR2, na.rm=TRUE)&n1[ ,4]<mean(BTR1, na.rm=TRUE),])
points(n2,pch=19,col=6,lwd=2,cex=3)

#Transition from Zone II to Zone IV
(n1=a[a[,1]>mean(BOR2, na.rm=TRUE)&a[,3]<mean(BOR1, na.rm=TRUE),])
(n2=n1[n1[ ,2]<mean(BTR2, na.rm=TRUE)&n1[ ,4]>mean(BTR1, na.rm=TRUE),])
points(n2,pch=1,lwd=2,cex=3)

#Transition from Zone III to Zone IV
(n1=a[a[,1]>mean(BOR2, na.rm=TRUE)&a[,3]>mean(BOR1, na.rm=TRUE),])
(n2=n1[n1[ ,2]<mean(BTR2, na.rm=TRUE)&n1[ ,4]>mean(BTR1, na.rm=TRUE),])
points(n2,pch=19,col="greenyellow",lwd=2,cex=3)

#Remainding in Zone IV
(n1=a[a[,1]>mean(BOR2, na.rm=TRUE)&a[,3]>mean(BOR1, na.rm=TRUE),])
(n2=n1[n1[ ,2]<mean(BTR2, na.rm=TRUE)&n1[ ,4]<mean(BTR1, na.rm=TRUE),])
points(n2,pch=1,lwd=2,cex=3)
text(BOR2,BTR2,1:length(BOR2))
}
