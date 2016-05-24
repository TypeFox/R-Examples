library("compositions")
data(SimulatedAmounts)
par(pch=20)
mydata <- simulateMissings(sa.groups5,dl=0.01,knownlimit=TRUE,MAR=0.05,MNARprob=0.05,SZprob=0.05)


cdata <- acomp(mydata)
cdata
mean(cdata)
mean(acomp(sa.groups5))
plot(acomp(sa.groups5))
plot(mean(cdata),add=T,col="red")
plot(mean(acomp(sa.groups5)),add=T,col="green")
mean(cdata - mean(cdata))
mm <-mean(cdata)
erg <-var(cdata)
print(erg)
svd(erg)
ellipses(mm,erg)
ellipses(mean(acomp(sa.groups5)),var(acomp(sa.groups5)))

cdata <- acomp(mydata)
plot(cdata)
plot(mean(cdata),add=T,col="blue")
plot(mean(acomp(sa.groups5)),add=T,col="green")
mean(cdata - mean(cdata))
mm <-mean(cdata)
#erg <-var(cdata)
svd(erg)
ellipses(mm,erg)
ellipses(mm,erg,r=2)



#cdata <- rcomp(mydata)
#cdata

#mean(cdata)
#mean(rcomp(sa.groups5))
#plot(rcomp(sa.groups5))
#plot(mean(cdata),add=T,col="red",pch=20)
#plot(mean(rcomp(sa.groups5)),add=T,col="green",pch=20)
#mean(cdata - mean(cdata)) # Nonsense because the difference is noncompositional


cdata <- aplus(mydata)
cdata
mean(cdata)
mean(aplus(sa.groups5))
plot(aplus(sa.groups5))
plot(mean(cdata),add=T,col="red",pch=20)
plot(mean(aplus(sa.groups5)),add=T,col="green",pch=20)
mean(cdata - mean(cdata))

cdata <- rplus(mydata)
cdata
mean(cdata)
mean(rplus(sa.groups5))
plot(rplus(sa.groups5))
plot(mean(cdata),add=T,col="red",pch=20)
plot(mean(rplus(sa.groups5)),add=T,col="green",pch=20)
mean(cdata - mean(cdata)) # Nonsense because the difference is non compositonal


plot(acomp(mydata))
plot(aplus(mydata))
plot(rcomp(mydata))
plot(rplus(mydata))

boxplot(acomp(mydata))
boxplot(aplus(mydata))
boxplot(rcomp(mydata))
boxplot(rplus(mydata))

barplot(acomp(mydata))
barplot(aplus(mydata))
barplot(rcomp(mydata))
barplot(rplus(mydata))
