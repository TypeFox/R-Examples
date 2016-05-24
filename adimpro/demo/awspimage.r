require(adimpro)
#
#     local polynomial smoothing of images
#
# construct an artificial example
isize <- 256
x<-seq(-1,1,length=isize)
fbi<-function(x,y,k,r){
z1<-sqrt(x^2+y^2)
theta<-asin(x/z1)
z<-sin(k*theta)
z[z1<r]<-sin(pi*z1/r)[z1<r]
z<-sign((x+y)*(x-y))*z
(z-min(z))/(max(z)-min(z))
}
xxx <- array(0,c(isize,isize,3))
xxx[,,1] <- .3+0.4*outer(x,x,"fbi",5,.5)
xxx[,,2] <- .25+0.5*outer(x,x,"fbi",7,.5)
xxx[,,3] <- .2+0.6*outer(x,x,"fbi",11,.5)
img0 <- make.image(xxx,gammatype="None")
show.image(img0,xaxt="n",yaxt="n")

xxx <- xxx+rnorm(xxx,0,.1)
xxx[xxx<0] <- 0
xxx[xxx>1] <- 1
img <- make.image(xxx,gammatype="None")
show.image(img,xaxt="n",yaxt="n")

X11(width=10,height=4)
imghat0 <- awsimage(img,hmax=12,graph=TRUE,varmodel="Constant")
X11(width=7,height=7)
imghat0b <- awsaniso(img,hmax=12,g=.1,rho=3,graph=TRUE,varmodel="Constant")
X11(width=10,height=4)
imghat1 <- awspimage(img,p=1,hmax=15,graph=TRUE,varmodel="Constant")
imghat2 <- awspimage(img,p=2,hmax=20,graph=TRUE,varmodel="Constant")

readline("Compare results:")

X11(width=9,height=6)
par(mfrow=c(2,3),mar=c(2,2,2,.2),mgp=c(2,1,0))
show.image(img0,xaxt="n",yaxt="n",main="True image")
show.image(img,xaxt="n",yaxt="n",main="Noisy image")
show.image(imghat0,xaxt="n",yaxt="n",main="Reconstruction (constant)")
show.image(imghat0b,xaxt="n",yaxt="n",main="Reconstruction (anisotrop constant)")
show.image(imghat1,xaxt="n",yaxt="n",main="Reconstruction (linear)")
show.image(imghat2,xaxt="n",yaxt="n",main="Reconstruction (quadratic)")

readline("Compare results (Edge enhancement):")

X11(width=9,height=6)
par(mfrow=c(2,3),mar=c(2,2,2,.2),mgp=c(2,1,0))
show.image(edges(img0,"Robertcross"),xaxt="n",yaxt="n",main="Edges true image")
show.image(edges(img,"Robertcross"),xaxt="n",yaxt="n",main="Edges noisy image")
show.image(edges(imghat0,"Robertcross"),xaxt="n",yaxt="n",main="Edges reconstruction (constant)")
show.image(edges(imghat0b,"Robertcross"),xaxt="n",yaxt="n",main="Edges reconstruction (anisotrop constant)")
show.image(edges(imghat1,"Robertcross"),xaxt="n",yaxt="n",main="Edges reconstruction (linear)")
show.image(edges(imghat2,"Robertcross"),xaxt="n",yaxt="n",main="Edges reconstruction (quadratic)")


z <- readline("End of demo, remove created objects (Y/N) :")

graphics.off()
if(toupper(z)!="N"){
rm(img,img0,imghat0,imghat0b,imghat1,imghat2, x, xxx, isize, fbi, z)
}
