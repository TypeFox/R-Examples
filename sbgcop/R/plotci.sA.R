"plotci.sA" <-
function(sA,ylabs=colnames(sA[,,1]),mgp=c(1.75,.75,0)) {
qA<-qM.sM(sA)
p<-dim(qA)[1]
tmp<-c(qA)
tmp<-tmp[tmp!=1]
par(mgp=mgp)

for(j in 1:p) {
plot(0,0,type="n",ylim=range(c(tmp),na.rm=TRUE),xlim=c(1,p),
     ylab=ylabs[j],xaxt="n", xlab="")
points( (1:p)[-j], qA[j,-j,2],pch=16,cex=.6 )
segments( x0=(1:p)[-j], y0=qA[j,-j,1], x1=(1:p)[-j], y1=qA[j,-j,3] )
abline(h=0,col="gray")
abline(v=j,col="gray")
               }
axis(side=1,at=1:p, labels=colnames(qA[,,1]),las=2)
                                              }

