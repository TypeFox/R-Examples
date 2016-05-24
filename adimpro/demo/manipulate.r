require(adimpro)
if(Sys.getenv("ImageMagick")==""){
cat("Please install ImageMagick\n")
} else {
#
#      changing the size of an image
#
X11(width=10.4,height=8)
par(mfrow=c(2,4),mar=c(1,1,1,.1),mgp=c(2,1,0))
img <- read.raw(system.file("img/rawimage.png",package="adimpro"))
readline("Choose a subimage by clipping:")

img1 <- clip.image(img,xaxt="n",yaxt="n")
readline("Shrink an image:")

img2 <- shrink.image(img,xt=200,yt=300)
show.image(img2,xaxt="n",yaxt="n",main="reduced size image")
summary(img2)
readline("Rotate an image:")

img3 <- rotate.image(img,angle=180)
show.image(img3,xaxt="n",yaxt="n",main="rotated image")
summary(img3)
readline("Change the color space of an image:")

img4 <- adjust.image(img,cspace="sRGB",gammatype="CIE",whitep="D55",
                     exposure=1.5)
show.image(img4,xaxt="n",yaxt="n",main="modified colors")
summary(img4)
readline("Create a gradient image:")

img5 <- imganiso2D(img)
show.image(img5,xaxt="n",yaxt="n",main="color coded image gradients")

readline("Create edge indicator images:")

img6 <- edges(img,abs=TRUE)
show.image(img6,xaxt="n",yaxt="n",main="Edge detection by abs(Laplacian)")

img7 <- edges(img,"Robertcross")
show.image(img7,xaxt="n",yaxt="n",main="Edge detection by Robertcross")

readline("Extract image data:")

ttt <- extract.image(img2)
cat("Range in red", range(ttt[,,1]), "Range in green", range(ttt[,,2]),
"Range in blue", range(ttt[,,3]),"\n")

readline("Extract image information:")

cat("Original image size",extract.info(img2,"Isize"),"\n")
z <- readline("End of demo, remove created objects (Y/N) :")

graphics.off()
if(toupper(z)!="N"){
rm(img, img1, img2, img3, img4, img5, img6, img7, ttt, z)
}
}


