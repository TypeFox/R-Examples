require(adimpro)
if(Sys.getenv("ImageMagick")==""){
cat("Please install ImageMagick\n")
} else {
system1 <- function(command){
# start a process in the background and return the process ID if .Platform$OS=="unix"
if(capabilities("X11")) {
system(paste(command,"&"))
if(.Platform$OS=="unix"){
# need a few milliseconds for the process to start
for(i in 1:1000) x <- runif(1)
# this's just to have some time before issuing the next command 
ttt <- strsplit(system(paste("ps aux|grep '",command,"'",sep=""),intern=TRUE)," +")
zzz <- grep("grep",ttt)
if(length(ttt[-zzz]) > 0) ttt[-zzz][[1]][2] else NULL
}
} else {
cat("Sorry no X-server on Windows\n")
""
}
}
#
#     I/O functions
#
#
#  Note: This demo requires ImageMagick
#
readline("Read an image in RAW-Format with conversion into RGB:")

# Images in camera RAW-format contain information as recorded by the 
# sensor of the digital camera. Usually only one out of three (four)
# color channels is available at each pixel. Which color is present 
# is coded in the Bayer mask. 
img <- read.raw(system.file("img/rawimage.png",package="adimpro"))
# Raw images can be read using function read.raw which uses dcraw  
# (Dave Coffin) for image conversion.
# Note: png images are used in package adimpro as a possibility to
# store RAW images in a non-camera specific format
X11(width=9,height=7)
par(mfrow=c(1,2),mar=c(2,2,2,.2),mgp=c(2,1,0))
show.image(img,xaxt="n",yaxt="n")
title("Part of an RAW image")
cat("Available information:\n")
summary(img)
#
#  Read and write other imaging formats
#
readline("Write image as png and read it again:")

write.image(img,"test.png",cspace="Adobe",gammatype="none")
#  Meta-Information on the image is stored as a comment.
#  This may be changed to using XMP-Profiles.
img2 <- read.image("test.png")
#  If Meta-Information is available in a comment it is used
#  to specify the properties of the adimpro-object
#
summary(img2)
show.image(img2,xaxt="n",yaxt="n")

#
#   write as 16-bit png
#
readline("Write image as png and reload image:")

write.image(img2,"testimage.png",gammatype="None",cspace="Adobe")
img3 <- read.image("testimage.png")
readline("Write image as png (cspace='Adobe', no gamma correction) and display result:")

write.image(img3,"testAdobe.png",gammatype="None",cspace="Adobe")
pid1 <- system1(paste(Sys.getenv("ImageMagick"),"display testAdobe.png",sep=""))
readline("Write image as png (cspace='Adobe', with gamma correction) and display result:")

write.image(img3,"testAdobegam.png",cspace="Adobe")
pid2 <- system1(paste(Sys.getenv("ImageMagick"),"display testAdobegam.png",sep=""))
readline("Write image as png (cspace='wGamut', no gamma correction) and display result:")

write.image(img3,"testwGamut.png",gammatype="None",cspace="wGamut")
pid3 <- system1(paste(Sys.getenv("ImageMagick"),"display testwGamut.png",sep=""))
readline("Write image as png (cspace='sRGB', no gamma correction) and display result:")

# This is lossy since sRGB has smaller gamut
write.image(img3,"testsRGB.png",gammatype="None",cspace="sRGB")
pid4 <- system1(paste(Sys.getenv("ImageMagick"),"display testsRGB.png",sep=""))
write.image(img3,"testsRGBgam8.jpg",cspace="sRGB",depth=8)
pid5 <- system1(paste(Sys.getenv("ImageMagick"),"display testsRGBgam8.jpg",sep=""))
readline("Write image as png (cspace='xyz') and display result:")
write.image(img3,"testxyz.png",cspace="xyz")

pid6 <- system1(paste(Sys.getenv("ImageMagick"),"display testxyz.png",sep=""))
readline("Write image as png (cspace='hsi') and display result:")

# This is lossy since HSI has smaller gamut
write.image(img3,"testhsi.png",cspace="hsi")
pid7 <- system1(paste(Sys.getenv("ImageMagick"),"display testhsi.png",sep=""))
readline("read these images and display results:")

if(.Platform$OS=="unix") system(paste("kill",pid1,pid2,pid3,pid4,pid5,pid6,pid7)) 
rm(pid1,pid2,pid3,pid4,pid5,pid6,pid7)
img4 <- read.image("testAdobe.png")
img5 <- read.image("testAdobegam.png")
img6 <- read.image("testwGamut.png")
img7 <- read.image("testsRGB.png")
img8 <- read.image("testxyz.png")
img9 <- read.image("testhsi.png")
img10 <- read.image("testsRGBgam8.jpg")
X11(width=9,height=7)
par(mfrow=c(2,4),mar=c(2,2,2,.2),mgp=c(2,1,0))
show.image(img2,xaxt="n",yaxt="n",main="Original")
show.image(img4,xaxt="n",yaxt="n",main="From testAdobe.png")
show.image(img5,xaxt="n",yaxt="n",main="From testAdobegam.png")
show.image(img6,xaxt="n",yaxt="n",main="From testwGamut.png")
show.image(img7,xaxt="n",yaxt="n",main="From testsRGB.png")
show.image(img8,xaxt="n",yaxt="n",main="From testxyz.png")
show.image(img9,xaxt="n",yaxt="n",main="From testhsi.png")
show.image(img10,xaxt="n",yaxt="n",main="From testsRGBgam8.jpg")
#
#  These images should look almost the same. Note that some degradation
#  occurs especially when gamma correction, a reduction to 8Bits,
#  a change of gamut and/or storage in a lossy format is involved.
#
readline("Show boxplots of image degradation due to gamma correction, a reduction to 8Bits,
  a change of gamut and/or storage in a lossy format:")

imgdata2 <- extract.image(img2)
imgdata4 <- extract.image(adjust.image(img4,cspace="Adobe",gammatype="None"))
imgdata7 <- extract.image(adjust.image(img7,cspace="Adobe",gammatype="None"))
imgdata10 <- extract.image(adjust.image(img10,cspace="Adobe",gammatype="None"))
X11(width=10,height=5)
par(mfrow=c(1,3),mar=c(2,2,2,.2),mgp=c(2,1,0))
boxplot(imgdata4-imgdata2)
title("Errors in testAdobe.png")
boxplot(imgdata7-imgdata2)
title("Errors in testsRGB.png (smaller Gamut)")
boxplot(imgdata10-imgdata2)
title("Errors in testsRGBgam8.jpg (smaller Gamut, 8bit, jpeg)")

z <- readline("End of demo, remove created objects (Y/N) :")

graphics.off()
if(toupper(z)!="N"){
for(i in c("test.png","testimage.png","testAdobe.png","testAdobegam.png","testwGamut.png","testsRGB.png","testxyz.png","testhsi.png","testsRGBgam8.jpg"))
file.remove(i)
rm(img,img2, img3, img4, img5, img6, img7, img8, img9, img10, imgdata2, imgdata4, imgdata7, imgdata10, z, system1)
}
}

