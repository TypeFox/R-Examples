library(compositions)
#library(compositions,lib.loc="../../compositions.Rcheck")
par(pch=20)
data(SimulatedAmounts)

geomSetA <- acomp(rbind(c(1,2,3),c(1,1,1)))
geomSetA
plot(geomSetA)
ellipses(acomp(geomSetA[1,]),ilrvar2clr(clrvar2ilr(diag(3))))

geomSetA <- acomp(rbind(c(a=1,b=2,c=3,d=4),c(1,1,1,1)))
geomSetA
plot(geomSetA,col=1:2,axes=list(side=-1:-3,pos=0.5))


ellipses(acomp(geomSetA[1,]),ilrvar2clr(clrvar2ilr(diag(4))))
ellipses(acomp(geomSetA[1,]),ilrvar2clr(clrvar2ilr(diag(4))),thinRatio=1)
ellipses(acomp(geomSetA[1,]),ilrvar2clr(clrvar2ilr(diag(4))),thinRatio=0)

plot(geomSetA,col=1:2,margin="b")
ellipses(acomp(geomSetA[1,]),ilrvar2clr(clrvar2ilr(diag(4))))

plot(rcomp(geomSetA),col=1:2,margin="b")
ellipses(rcomp(geomSetA[1,]),ilrvar2clr(clrvar2ilr(diag(4))))


geomSetA <- acomp(rbind( c(2,1,1),c(1,2,1),c(1,1,2),c(1,1,1)))
delta <- acomp(c(0.4,0.9,3))
plot(geomSetA+delta)
straight(acomp(c(2,1,1))+delta,acomp(c(1,2,1))-acomp(c(2,1,1)))
ellipses(acomp(c(1,1,2))+delta,var=diag(3),r=norm(acomp(c(1,1,2))-acomp(c(1,1,1))))
ellipses(mean(geomSetA+delta),var=var(geomSetA+delta),r=2)


plot(rcomp(geomSetA+delta))
straight(rcomp(acomp(c(2,1,1))+delta),rcomp(acomp(c(1,2,1))+delta)-rcomp(acomp(c(2,1,1))+delta))
ellipses(rcomp(acomp(c(1,1,2)+delta)),var=diag(3)-1/3,r=norm(rcomp(acomp(c(1,1,2))+delta)-rcomp(acomp(c(2,1,1))+delta)))
ellipses(mean(rcomp(geomSetA+delta)),var=var(rcomp(geomSetA+delta)),r=2)

geomSetA <- acomp(rbind( c(2,1,1,1),c(1,2,1,1),c(1,1,2,1),c(1,1,1,2),c(1,1,1,1)))
delta <- acomp(c(0.4,0.9,3,1.2))
#delta <- acomp(c(1,1,1,1))
plot(geomSetA+delta,col=1:5,pch=1:5)
straight(acomp(c(2,1,1,1))+delta,acomp(c(1,2,1,1))-acomp(c(2,1,1,1)),steps=100)
#straight(acomp(c(2,1,1,1))+delta,acomp(c(1,2,1,1))-acomp(c(2,1,1,1)),steps=10)
#ellipses(acomp(c(1,1,2,1))+delta,var=diag(4)-1/4,r=norm(acomp(c(1,1,2,1))-acomp(c(1,1,1,1))))
replot(onlyPanel=c(3,1))
segments.acomp(acomp(c(2,1,1,1))+delta,acomp(c(1,2,1,1))+delta,col="red")

ellipses(acomp(c(1,1,1,1))+delta,var=diag(4)-1/4,r=norm(acomp(c(1,1,2,1))-acomp(c(1,1,1,1))),col="blue")
ellipses(mean(geomSetA+delta),var=var(geomSetA+delta),r=2,col="red")

ellipses(mean(geomSetA+delta),var=var(geomSetA+delta),r=2,col="red",thinRatio=0.1)

plot(rcomp(geomSetA+delta),col=1:5,pch=1:5)
straight(rcomp(acomp(c(2,1,1,1))+delta),rcomp(acomp(c(1,2,1,1))+delta)-rcomp(acomp(c(2,1,1,1))+delta))
ellipses(rcomp(acomp(c(1,1,2,1)+delta)),var=diag(4)-1/4,r=norm(rcomp(acomp(c(1,1,2,1))+delta)-rcomp(acomp(c(2,1,1,1))+delta)))
ellipses(mean(rcomp(geomSetA+delta)),var=var(rcomp(geomSetA+delta)),r=2)

geomSetA <- aplus(rbind( c(2,1,1),c(1,2,1),c(1,1,2),c(1,1,1),c(1.5,1.5,1.5)))
delta <- aplus(c(0.4,0.9,3))
#delta <- aplus(c(1,1,1))
plot(geomSetA+delta,xlim=c(0.1,10),ylim=c(0.1,10),logscale=TRUE)
straight(aplus(c(2,1,1))+delta,aplus(c(1,2,1))-aplus(c(2,1,1)))
ellipses(aplus(c(1.5,1.5,1.5))+delta,var=diag(3),r=norm(aplus(c(1,1,2))-aplus(c(1,1,1))))
ellipses(aplus(c(1,1,1))+delta,var=diag(3),r=norm(aplus(c(1,1,2))-aplus(c(1,1,1))))
ellipses(mean(geomSetA+delta),var=var(geomSetA+delta),r=2)
lines(geomSetA+delta)
geomSetA <- rplus(rbind( c(2,1,1),c(1,2,1),c(1,1,2),c(1,1,1),c(1.5,1.5,1.5)))
delta <- rplus(c(0.4,0.9,3))
plot(geomSetA+delta,xlim=c(0.1,10),ylim=c(0.1,10))
#straight(aplus(c(2,1,1))+delta,aplus(c(1,2,1))-aplus(c(2,1,1)))
ellipses(rplus(c(1.5,1.5,1.5))+delta,var=diag(3),r=norm(rplus(c(1,1,2))-rplus(c(1,1,1))))
ellipses(mean(geomSetA+delta),var=var(geomSetA+delta),r=2)

plot(rcomp(sa.lognormals))
ternaryAxis(side=1:3,pos=0,col.axis="red",col.lab="green")
ternaryAxis(side=1:3,at=1:9/10,labels=expression(9:1,4:1,7:3,3:2,1:1,2:3,3:7,1:4,1:9),pos=0,col.axis="red",col.lab="green")
ternaryAxis(side=rep(-1:-3,3),labels=paste(seq(20,80,by=20),"%"),
            pos=rep(c(0,0.5,1),each=3),col.axis=1:3,col.lab="green")
ternaryAxis(side=rep(1:3,3),at=1:9/10,labels=expression(9:1,4:1,7:3,3:2,1:1,2:3,3:7,1:4,1:9),pos=rep(c(0,0.5,1),each=3))
