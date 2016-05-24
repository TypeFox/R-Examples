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
cat("Segmentation with specified value of level\n")
simg1 <- segment(img,level=14600,delta=3000,hmax=20,select=FALSE,graph=TRUE)
cat("Segmentation with selection of a homogeneous region and 
    extraction of a connected set\n")
simg2 <- segment(img,hmax=20,select=TRUE,connect=TRUE)
show.image(simg2)
cat("Extract connected set to create a mask\n")
simg2data <- extract.image(simg2)
levels <- as.integer(names(table(simg2data)))
mask <- simg2data==levels[length(levels)]
image(mask)
z <- readline("End of demo, remove created objects (Y/N) :")

graphics.off()
if(toupper(z)!="N"){
rm(img,simg1,simg2,simg2data,levels,mask,z)
}
}
