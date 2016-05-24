### R code from vignette source 'alphashape3d.rnw'

###################################################
### code chunk number 1: alphashape3d.rnw:102-108
###################################################
library(alphahull)
library(rgl)
library(geometry)
#library(penalizedSVM)
library(alphashape3d)
options(prompt = "R> ", continue = "+  ")


###################################################
### code chunk number 2: alphashape3d.rnw:299-328
###################################################
n<-300
m<-0
data<-matrix(0,n,2)
while(m<n){
x<-runif(1)
y<-runif(1)
if(((x-0.5)^2+(y-0.5)^2<=(0.5)^2 )&((x-0.5)^2+(y-0.5)^2>=(0.25)^2)){
m<-m+1
data[m,]<-c(x,y)
}
}
x<-data
par(mfrow=c(1,3))
alpha1=0.04
alpha2=0.2
alpha3=1
plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main=expression(paste(alpha," = 0.04 ")),axes=F,cex.main=2.7,asp=TRUE)
axis(1,cex.axis=2.7)
axis(2,cex.axis=2.7)
plot(ahull(x,alpha=alpha1),do.sha=TRUE,wpoints=TRUE,col=c(6,4,rep(1,3)),lwd=2,add=T)
plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main=expression(paste(alpha," = 0.2 ")),axes=F,cex.main=2.7,asp=TRUE)
axis(1,cex.axis=2.7)
axis(2,cex.axis=2.7)
plot(ahull(x,alpha=alpha2),do.sha=TRUE,wpoints=TRUE,col=c(6,4,rep(1,3)),lwd=2,add=T)
plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main=expression(paste(alpha," = 1 ")),axes=F,cex.main=2.7,asp=TRUE)
axis(1,cex.axis=2.7)
axis(2,cex.axis=2.7)
plot(ahull(x,alpha=alpha3),do.sha=TRUE,wpoints=TRUE,col=c(6,4,rep(1,3)),lwd=2,add=T)
par(mfrow=c(1,1))


###################################################
### code chunk number 3: alphashape3d.rnw:364-384
###################################################
n<-800
n2<-20
t<-runif(n,0,pi/1.2)
x<-sin(6*t)
y<-cos(6*t)
z<-t
x1<-matrix(0,nr=n,ncol=n2)
y1<-matrix(0,nr=n,ncol=n2)
z1<-matrix(0,nr=n,ncol=n2)
for(i in 1:n){
norm<-matrix(rnorm(n2*3),nc=3)
aux<-norm/sqrt(rowSums(norm^2))
r<-runif(n2,0,0.15)
aux<-r*aux
x1[i,]<-aux[,1]+x[i]
y1[i,]<-aux[,2]+y[i]
z1[i,]<-aux[,3]+z[i]
}
x<-cbind(as.numeric(t(x1)),as.numeric(t(y1)),as.numeric(t(z1)))
x<-rotate3d(x,pi/2.8,0.9,0.9,0.5)


###################################################
### code chunk number 4: alphashape3d.rnw:388-391
###################################################
bg3d("white")
as<-ashape3d(x,alpha=0.03)
plot(as,col=c(4,4,4))


###################################################
### code chunk number 5: alphashape3d.rnw:393-396
###################################################
as<-ashape3d(x,alpha=0.35)
bg3d("white")
plot.ashape3d(as,col=c(4,4,4))


###################################################
### code chunk number 6: alphashape3d.rnw:398-401
###################################################
as<-ashape3d(x,alpha=0.5)
bg3d("white")
plot.ashape3d(as,col=c(4,4,4))


###################################################
### code chunk number 7: alphashape3d.rnw:469-472
###################################################
x1<-c(0.5915,0.6230,0.9689,0.8248,0.9392,0.8156,0.2050,0.9757,0.0957,0.4139)
y1<-c(0.472,0.619,0.304,0.197,0.716,0.575,0.507,0.574,0.996,0.893)
x<-cbind(x1,y1)


###################################################
### code chunk number 8: alphashape3d.rnw:478-481
###################################################
par(mfrow=c(1,1))
print(plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F))
print(plot(delvor(x),col=1:3,xlab="",ylab="",add=T))


###################################################
### code chunk number 9: alphashape3d.rnw:483-532
###################################################
# Delaunay triangulation in R3 and circumsphere
# This function is to choose a not too large circumsphere
# We consider the circumspheres contained in [0,1]^3
DTcube<-function(v){
return(apply(cbind(v[,1],1-v[,1],v[,2],1-v[,2],v[,3],1-v[,3]),1,min))
}
# x..........data
n<-50
x<-cbind(runif(n),runif(n),runif(n))
# delaunayn return a m*4 matrix (m tetrahedra with 4 indexed vertices abcd)
tc<-delaunayn(x)
# Coordenates a
v1<-x[tc[,1],]
v2<-x[tc[,2],]
v3<-x[tc[,3],]
v4<-x[tc[,4],]

nv1<-apply(v1^2,1,sum)
nv2<-apply(v2^2,1,sum)
nv3<-apply(v3^2,1,sum)
nv4<-apply(v4^2,1,sum)


Dx<-numeric()
Dy<-numeric()
Dz<-numeric()
ct<-numeric()
a<-numeric()
for (i in 1:dim(v1)[1]){
tetra<-rbind(v1[i,],v2[i,],v3[i,],v4[i,])
nor<-c(nv1[i],nv2[i],nv3[i],nv4[i])
Dx[i]<-det(cbind(nor,tetra[,2:3],rep(1,4)))
Dy[i]<--det(cbind(nor,tetra[,c(1,3)],rep(1,4)))
Dz[i]<-det(cbind(nor,tetra[,1:2],rep(1,4)))
ct[i]<-det(cbind(nor,tetra[,1:3]))
a[i]<-det(cbind(tetra[,1:3],rep(1,4)))
}
# Circumspheres with centers (ctx,cty,ctz) and radius rad
ctx<-0.5*Dx/a
cty<-0.5*Dy/a
ctz<-0.5*Dz/a
rad<-0.5*sqrt(Dx^2+Dy^2+Dz^2-4*a*ct)/abs(a)
# We determine which circumspheres are contained in [0,1]^3
d1<-DTcube(cbind(ctx,cty,ctz))
inS<-rad<d1
# we choose an "intermediate" one to draw (not too large, not too big)
# cualSp1<-which(rad==max(rad[inS]))
cualSp1<-which(rad==sort(rad[inS])[floor(length(rad[inS])/2)])
sp1<-c(ctx[cualSp1],cty[cualSp1],ctz[cualSp1],rad[cualSp1])


###################################################
### code chunk number 10: alphashape3d.rnw:534-541
###################################################
rgl.viewpoint(60)
rgl.light(120,60)
rgl.bg(col="white")
tetramesh(tc,x,alpha=0.1,color="red")  
tetramesh(matrix(tc[cualSp1,],nr=1,nc=4),x, alpha=1,color="red",add=TRUE,clear=FALSE)
rgl.spheres(sp1[1],sp1[2],sp1[3],sp1[4],alpha=0.5,color="gray")
points3d(x[,1],x[,2],x[,3],color=1)


###################################################
### code chunk number 11: alphashape3d.rnw:578-583
###################################################
n<- 4000
r1<- 0.5
r2<- 2
alp1<- 0.3
alp2<- 0.5


###################################################
### code chunk number 12: alphashape3d.rnw:590-595
###################################################
n <- 4000
T1 <- rtorus(n/2, 0.5, 2)
T2 <- rtorus(n/2, 0.5, 2, ct = c(2, 0, 0), rotx = pi/3)
x <- rbind(T1, T2)
points3d(x, col = 4)


###################################################
### code chunk number 13: alphashape3d.rnw:600-604
###################################################
rgl.viewpoint(20)
rgl.light(120,60)
bg3d("white")
rgl.points(x, col = 4)


###################################################
### code chunk number 14: alphashape3d.rnw:628-633 (eval = FALSE)
###################################################
## n <- 4000
## T1 <- rtorus(n/2, 0.5, 2)
## T2 <- rtorus(n/2, 0.5, 2, ct = c(2, 0, 0), rotx = pi/3)
## x <- rbind(T1, T2)
## points3d(x, col = 4)


###################################################
### code chunk number 15: alphashape3d.rnw:635-636
###################################################
rgl.close()


###################################################
### code chunk number 16: alphashape3d.rnw:641-644
###################################################
alphashape3d <- ashape3d(x, alpha = c(0.3, 0.5))
class(alphashape3d)
names(alphashape3d)


###################################################
### code chunk number 17: alphashape3d.rnw:656-657
###################################################
head(alphashape3d$tetra)


###################################################
### code chunk number 18: alphashape3d.rnw:664-665
###################################################
head(alphashape3d$triang)


###################################################
### code chunk number 19: alphashape3d.rnw:670-671
###################################################
head(alphashape3d$edge)


###################################################
### code chunk number 20: alphashape3d.rnw:678-679
###################################################
head(alphashape3d$vertex)


###################################################
### code chunk number 21: alphashape3d.rnw:692-696
###################################################
rgl.viewpoint(20)
rgl.light(120, 60)
bg3d("white")
plot(alphashape3d, col = c(4, 4, 4), indexAlpha = "all")


###################################################
### code chunk number 22: alphashape3d.rnw:698-699
###################################################
rgl.close()


###################################################
### code chunk number 23: alphashape3d.rnw:712-716
###################################################
rgl.viewpoint(20)
rgl.light(120,60)
bg3d("white")
plot.ashape3d(alphashape3d,indexAlpha=1, col = rep(4,3))


###################################################
### code chunk number 24: alphashape3d.rnw:719-723
###################################################
rgl.viewpoint(20)
rgl.light(120,60)
bg3d("white")
plot.ashape3d(alphashape3d,indexAlpha=2, col = rep(4,3))


###################################################
### code chunk number 25: alphashape3d.rnw:765-766
###################################################
comp <- components_ashape3d(alphashape3d, indexAlpha = "all")


###################################################
### code chunk number 26: alphashape3d.rnw:771-773
###################################################
table(comp[[1]])
table(comp[[2]])


###################################################
### code chunk number 27: alphashape3d.rnw:781-785
###################################################
rgl.viewpoint(20)
rgl.light(120, 60)
bg3d("white")
plot(alphashape3d, byComponents = TRUE, indexAlpha = "all")


###################################################
### code chunk number 28: alphashape3d.rnw:787-788
###################################################
rgl.close()


###################################################
### code chunk number 29: alphashape3d.rnw:795-799
###################################################
rgl.viewpoint(20)
rgl.light(120,60)
bg3d("white")
plot(alphashape3d, byComponents = TRUE, indexAlpha = 1)


###################################################
### code chunk number 30: alphashape3d.rnw:801-805
###################################################
rgl.viewpoint(20)
rgl.light(120,60)
bg3d("white")
plot(alphashape3d, byComponents = TRUE, indexAlpha = 2)


###################################################
### code chunk number 31: alphashape3d.rnw:819-820
###################################################
volume_ashape3d(alphashape3d, indexAlpha = 1)


###################################################
### code chunk number 32: alphashape3d.rnw:822-823
###################################################
volt<-volume_ashape3d(alphashape3d, indexAlpha = 1)


###################################################
### code chunk number 33: alphashape3d.rnw:828-829
###################################################
volume_ashape3d(alphashape3d, indexAlpha = 1, byComponents = TRUE)


###################################################
### code chunk number 34: alphashape3d.rnw:840-841
###################################################
inashape3d(alphashape3d, indexAlpha = "all", c(-2.5, 2.5, 0))


###################################################
### code chunk number 35: alphashape3d.rnw:846-855
###################################################
np <- 50
x <- seq(-2.5, 2.5, length = np)
points <- cbind(x, -x, rep(0, np))
in3d <- inashape3d(alphashape3d, indexAlpha = 1, points = points)
rgl.viewpoint(20)
bg3d("white")
plot(alphashape3d, col = c(4, 4, 4), indexAlpha = 1, transparency = 0.2)
colours <- ifelse(in3d, "green", "red")
rgl.points(points, col = colours)


###################################################
### code chunk number 36: alphashape3d.rnw:867-872
###################################################
rgl.viewpoint(20)
bg3d("white")
plot(alphashape3d, indexAlpha = 1, col = c(4, 4, 4), transparency = 0.2)
colours <- ifelse(in3d, "green", "red")
rgl.points(points, col = colours)


###################################################
### code chunk number 37: alphashape3d.rnw:903-908
###################################################
normv <- surfaceNormals(alphashape3d, indexAlpha = 1)
class(normv)
names(normv)
head(normv$normals)
head(normv$centers) 


###################################################
### code chunk number 38: alphashape3d.rnw:928-938
###################################################
offFilename <- "alphashape3d.off"
write("OFF", file = offFilename, append = FALSE, sep = " ")
aux <- alphashape3d$triang
tr <- aux[aux[, 9] == 2 | aux[, 9]==3, c("tr1", "tr2", "tr3")]
write(c(dim(alphashape3d$x)[1], dim(tr)[1], 0), file = offFilename, 
append = TRUE, sep = " ")
write(t(alphashape3d$x), file = offFilename, ncolumns = 3, 
append = TRUE, sep = " ")
write(t(cbind(rep(3, dim(tr)[1]), tr - 1)), 
file = offFilename, ncolumns = 4, append = TRUE, sep = " ")


