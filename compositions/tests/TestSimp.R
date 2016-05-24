library("compositions")

#js <- read.table("juraset.dat",skip=17,header=TRUE)
#js$Land <- factor( c("Wald","Weide","Wiese","Acker")[js$Land] )
#js$Rock <- factor( c("Argovian","Kimmeridgian","Sequanian","Portlandian","Quaternary")[js$Rock] )

#Rock <- js$Rock
#Land <- js$Land

#cdata <- js[,c("Cd","Cu","Pb","Co","Cr","Ni","Zn")]
#cdata <- data.matrix(js[,c("Cd","Cu","Pb","Co","Cr","Ni","Zn")])

#cd1 <- cdata[,1:3]
#cd2 <- cdata[,4:7]
data(SimulatedAmounts)
cdata <- sa.groups5
Land  <- sa.groups5.area
cd1 <- cdata[,1:3]
cd2 <- cdata[,4:5]

# Transformations
# clo 
clo(c(1,2,3))
clo(matrix(1:4,ncol=2))
data(iris)
clo(iris[,1:4])
clo(0.5)
clo(matrix(0.5))
clo(matrix(0.5,nrow=5))
clo(matrix(0.5,ncol=5))

clo(iris,c("Sepal.Length","Petal.Length"))
clo(iris,c(2,3))

checker <- function(x,y) {
  x<-unclass(x)
  y<-unclass(y)
  if( sum(c(x-y)^2) > 1E-10 )
    stop("Wrong results in ", deparse(substitute(x)))
  x
}
clr(cdata)
checker( clrInv(clr(cdata)) , clo(cdata) )
ilr(cdata)
checker( ilrInv(ilr(cdata)) , clo(cdata) )
alr(cdata)
checker( alrInv(alr(cdata)) , clo(cdata) )
cpt(cdata)
checker( cptInv(cpt(cdata)) , clo(cdata) )
ipt(cdata)
checker( iptInv(ipt(cdata)) , clo(cdata) )
apt(cdata)
checker( aptInv(apt(cdata)) , clo(cdata) )
ilt(cdata)
checker( iltInv(ilt(cdata)) , cdata )
iit(cdata)
checker( iitInv(iit(cdata)) , cdata )

clr(c(a=1,2,3))
ilr(c(a=1,2,3))
alr(c(a=1,2,3))
cpt(c(a=1,2,3))
ipt(c(a=1,2,3))
apt(c(a=1,2,3))
ilt(c(a=1,2,3))
iit(c(a=1,2,3))
checker( clrInv(clr(c(a=1,2,3))) , clo(c(1,2,3)))
checker( ilrInv(ilr(c(a=1,2,3))) , clo(c(1,2,3)))
checker( alrInv(alr(c(a=1,2,3))) , clo(c(1,2,3)))
checker( cptInv(cpt(c(a=1,2,3))) , clo(c(1,2,3)))
checker( iptInv(ipt(c(a=1,2,3))) , clo(c(1,2,3)))
checker( aptInv(apt(c(a=1,2,3))) , clo(c(1,2,3)))
checker( iltInv(ilt(c(a=1,2,3))) , c(1,2,3))
checker( iitInv(iit(c(a=1,2,3))) , c(1,2,3))

# mean

mean(acomp(cdata))
mean(rcomp(cdata))
mean(aplus(cdata))
mean(rplus(cdata))

meanCol(cdata)
meanCol(clo(cdata))
clo(meanCol(cdata))

# var (Variation Matrix)
var(rcomp(cdata))
var(acomp(cdata))
var(aplus(cdata))
var(rplus(cdata))

# clr
clr(mean(acomp(cdata)))
meanCol(clr(cdata))
clrInv(meanCol(clr(cdata)))
mean(acomp(cdata))

# ilr
ilr(mean(acomp(cdata)))
meanCol(ilr(cdata))
ilrInv(meanCol(ilr(cdata)))
mean(acomp(cdata))

# alr
alr(mean(acomp(cdata)))
meanCol(alr(cdata))
alrInv(meanCol(alr(cdata)))
mean(acomp(cdata))

# Operations
mean(acomp(3 * (cdata - mean(acomp(cdata)))))


# barplot
barplot(acomp(cdata[1:10,]))
barplot(rcomp(cdata[1:10,]))
barplot(aplus(cdata[1:10,]))
barplot(rplus(cdata[1:10,]))

barplot(mean(acomp(cdata)))
barplot(mean(rcomp(cdata)))
barplot(mean(aplus(cdata)))
barplot(mean(rplus(cdata)))

# piechart
pie(mean(acomp(cdata)))
pie(mean(rcomp(cdata)))
pie(mean(aplus(cdata)))



# Triangular Diagrams
plot(acomp(cdata[,1:3]))
mean(acomp(cdata[,1:3]))
plot(acomp(cdata),margin="rcomp",pca=TRUE)
plot(acomp(cdata),margin="acomp",pca=TRUE)
#plot(acomp(cdata),margin="Cd",pca=TRUE) # Bug!!!
mean(acomp(cdata))

boxplot(acomp(cdata))              # boxplotscale
boxplot(acomp(cdata),Land,notch=TRUE) # notch
#boxplot(acomp(cdata),Rock)

boxplot(acomp(cdata),log=FALSE)
boxplot(acomp(cdata),Land,log=FALSE)
#boxplot(acomp(cdata),Rock,log=FALSE)

boxplot(acomp(cdata),log=FALSE,ylim=c(0,5))
boxplot(acomp(cdata),Land,log=FALSE,ylim=c(0,5))
#boxplot(acomp(cdata),Rock,log=FALSE,ylim=c(0,5))

qqnorm(acomp(cdata),alpha=100)
qqnorm(acomp(cdata),alpha=0.05)
qqnorm(acomp(cdata[,-3]),alpha=0.05)
qqnorm(acomp(cdata[,-3]),alpha=0.05)
 
#boxplot(acomp(cdata[,-1]),js$Cd)
plot(Land,data.matrix(cdata)%*% rep(1,ncol(cdata)))


# rcomp.plots

boxplot(rcomp(cdata))
boxplot(rcomp(cdata),Land)
#boxplot(rcomp)(cdata,Rock)

qqnorm(rcomp(cdata))
qqnorm(rcomp(cdata),alpha=0.05)
qqnorm(rcomp(cdata[,1:3]),alpha=0.05)
plot(acomp(cdata[,1:3]))
ellipses(mean(acomp(cdata[,1:3])), var(acomp(cdata,1:3)),col="red",r=2)

ellipses(mean(rcomp(cdata[,1:3])), var(rcomp(cdata[,1:3])),col="blue",r=2)


plot(rplus(cdata[,1:2]))
ellipses(rplus(mean(rplus(cdata[,1:2]))), var(rplus(cdata[,1:2])),col="blue",r=2)
ellipses(aplus(mean(aplus(cdata[,1:2]))), var(aplus(cdata[,1:2])),col="red",r=2)

plot(aplus(cdata[,1:2]))
ellipses(aplus(mean(aplus(cdata[,1:2]))), var(aplus(cdata[,1:2])),col="red",r=2)
ellipses(rplus(mean(rplus(cdata[,1:2]))), var(rplus(cdata[,1:2])),col="blue",r=2)


straight(acomp(c(1,1,1)),c(2,1,3))


boxplot(rcomp(cdata[,-1]),cdata[,"Cd"])

# biplot, princomp

biplot(princomp(cdata))
biplot(princomp(acomp(cdata)))
biplot(princomp(rcomp(cdata)))
biplot(princomp(aplus(cdata)),choice=c(2,3))
biplot(princomp(rplus(cdata)),choice=c(2,3))

summary(princomp(cdata))
summary(princomp(acomp(cdata)))
summary(princomp(rcomp(cdata)))
summary(princomp(aplus(cdata)))
summary(princomp(rplus(cdata)))

loadings(princomp(cdata))
loadings(princomp(acomp(cdata)))
loadings(princomp(rcomp(cdata)))
loadings(princomp(aplus(cdata)))
loadings(princomp(rplus(cdata)))


#names
meanCol(cdata[,1:3])
mean(acomp(cdata[,1:3]))
oneOrDataset(c(a=1,b=2,c=3))

# covariance
cov(acomp(cd1),acomp(cd2))
cov(rcomp(cd1),rcomp(cd2))
cov(aplus(cd1),aplus(cd2))
cov(rplus(cd1),rplus(cd2))

#tmp <- princov(acomp(cd1,cd2))
#tmp
#plot(tmp)
#biplot(tmp)


tmp <- princomp(acomp(cdata))
tmp
plot(tmp)
biplot(tmp)   # Fehler !!!
biplot(tmp,ch=2:3)

#pplot(acomp(cdata))

# Simulated data



# Loading data
  read.geoeas("TRUE.DAT")
  read.geoEAS("TRUE.DAT")
 # read.standard("LAKE.txt")
