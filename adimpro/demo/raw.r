require(adimpro)
if(Sys.getenv("ImageMagick")==""){
cat("Please install ImageMagick\n")
} else {
#
#     I/O functions
#
#
#  Note: This demo requires ImageMagick
#
readline("Read an image in RAW-Format:")

# Images in camera RAW-format contain information as recorded by the 
# sensor of the digital camera. Usually only one out of three (four)
# color channels is available at each pixel. Which color is present 
# is coded in the Bayer mask. 
img <- read.raw(system.file("img/rawimage.png",package="adimpro"),type="RAW")
# Raw images can be read using function read.raw which uses dcraw  
# (Dave Coffin) for image conversion.
# Note: png images are used in package adimpro as a possibility to
# store RAW images in a non-camera specific format
X11(width=10,height=7)
par(mfrow=c(2,3),mar=c(2,2,2,.2),mgp=c(2,1,0))
show.image(img,xaxt="n",yaxt="n")
title("Part of an RAW image")
cat("Available information:\n")
summary(img)
readline("Convert image to RGB:")

img1 <- develop.raw(img,method="HALF")
show.image(img1,xaxt="n",yaxt="n")
title("Half-size color image")
img2 <- develop.raw(img,method="BILINEAR")
show.image(img2,xaxt="n",yaxt="n")
title("Bilinear interpolation")
readline("Clipping an image in RAW-Format:")

img3 <- clip.image(img)
readline("Writing an image in RAW-Format:")

write.raw(img3,"rawtest")
# Note: This is how img/rawimage.png was generated from the RAW image
# CRW_2844.CRW (Canon EOS 300D DIGITAL). Information available with 
# the original image are stored as a comment in the generated png.

readline("Read an image in RAW-Format with conversion into RGB:")

img4 <- read.raw("rawtest.png")
show.image(img4,xaxt="n",yaxt="n")
summary(img4)

readline("Smooth and demosaic the RAW image:")
img5 <- awsraw(img,hmax=10,graph=TRUE)
X11(width=12,height=5)
par(mfrow=c(1,3),mar=c(1,1,3,.5),mgp=c(2,1,0))
show.image(img,xaxt="n",yaxt="n")
title("RAW image")
show.image(img1,xaxt="n",yaxt="n")
title("Demosaicing only")
show.image(img5,xaxt="n",yaxt="n")
title("Smoothing and demosaicing")
summary(img5)

z <- readline("End of demo, remove created objects (Y/N) :")

graphics.off()
if(toupper(z)!="N"){
file.remove("rawtest.png")
rm(img,img1,img2,img3,img4,z)
}
}