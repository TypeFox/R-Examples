#2D example : female and male Gorillas (cf. Dryden and Mardia, 1998)

data(gorf.dat)
data(gorm.dat)

n1<-dim(gorf.dat)[3]
n2<-dim(gorm.dat)[3]
k<-dim(gorf.dat)[1]
m<-dim(gorf.dat)[2]
gor.dat<-array(0,c(k,2,n1+n2))
gor.dat[,,1:n1]<-gorf.dat
gor.dat[,,(n1+1):(n1+n2)]<-gorm.dat

plotshapes(gorf.dat,gorm.dat)
gorf<-procGPA(gorf.dat)
gorm<-procGPA(gorm.dat)
plotshapes(gorf$rotated,gorm$rotated)

gor<-procGPA(gor.dat)
shapepca(gor,type="r",mag=3)
cat("First 3 PCs for pooled gorilla data")

for (ii in 1:1000000){a<-1} 

shapepca(gor,type="v",mag=3)
cat("First 3 PCs for pooled gorilla data")

for (ii in 1:1000000){a<-1} 



gor.gp<-c(rep("f",times=30),rep("m",times=29))
x<-cbind(gor$size,gor$rho,gor$scores[,1:3])
pairs(x,panel=function(x,y) text(x,y,gor.gp),
   label=c("s","rho","score 1","score 2","score 3"))


#Some tests                 
test1<-testmeanshapes(gorf.dat,gorm.dat)


print(test1)


cat("Highly significant difference in mean shape")

gorf<-procGPA(gorf.dat)
gorm<-procGPA(gorm.dat)
mag<-2
TT<-gorf$mshape
YY<-gorm$mshape
par(mfrow=c(1,2))
YY<-TT+(YY-TT)*mag
tpsgrid(TT,YY,-150,-150,300,2,0.1,22) 
title("TPS grid: Female mean (left) to Male mean (right)")   

cat("##########################################################\n")
cat("#3D example\n#")
cat("##########################################################\n")

# Male macaques data (Dryden/Mardia 1998) 

data(macm.dat)
out<-procGPA(macm.dat)

par(mfrow=c(2,2))
plot(out$rawscores[,1],out$rawscores[,2],xlab="PC1",ylab="PC2")
title("PC scores")
plot(out$rawscores[,2],out$rawscores[,3],xlab="PC2",ylab="PC3")
plot(out$rawscores[,1],out$rawscores[,3],xlab="PC1",ylab="PC3")
plot(out$size,out$rho,xlab="size",ylab="rho")
title("Size versus shape distance")

