require(adimpro)
if(Sys.getenv("ImageMagick")==""){
cat("Please install ImageMagick\n")
} else {
#
#     smoothing images
#
X11(width=8,height=4)
par(mfrow=c(1,2),mar=c(2,2,2,.2),mgp=c(2,1,0))
img <- read.image(system.file("img/wias-noise.ppm",package="adimpro"))
img0 <- read.image(system.file("img/wias.ppm",package="adimpro"))
show.image(img0,main="True image",xaxt="n",yaxt="n")
show.image(img,main="Noisy image",xaxt="n",yaxt="n")
readline("Structural adaptive isotropic smoothing:")

X11(width=10,height=4)
imghat <- awsimage(img,hmax=12,graph=TRUE)
readline("Structural adaptive anisotropic smoothing:")

X11(width=7,height=7)
imghat2 <- awsaniso(img,hmax=12,graph=TRUE,g=.1,rho=3)
readline("Compare results:")

X11(width=7,height=7)
par(mfrow=c(2,2),mar=c(2,2,2,.2))
show.image(img,xaxt="n",yaxt="n")
title("Noisy Original")
show.image(img0,xaxt="n",yaxt="n")
title("True image")
show.image(imghat,xaxt="n",yaxt="n")
title("Isotropic adaptive")
show.image(imghat2,xaxt="n",yaxt="n")
title("Anisotropic adaptive")
readline("Edge detection:")

X11(width=7,height=7)
par(mfrow=c(2,2),mar=c(2,2,2,.2))
show.image(edges(img,abs=TRUE),xaxt="n",yaxt="n")
title("Edges Noisy Original")
show.image(edges(imghat,abs=TRUE),xaxt="n",yaxt="n")
title("Edges Isotropic adaptive")
show.image(edges(imghat2,abs=TRUE),xaxt="n",yaxt="n")
title("Edges Anisotropic adaptive")
show.image(edges(img0,abs=TRUE),xaxt="n",yaxt="n")
title("Edges true image")


z <- readline("End of demo, remove created objects (Y/N) :")

graphics.off()
if(toupper(z)!="N"){
rm(img,img0,imghat,imghat2, z)
}
}