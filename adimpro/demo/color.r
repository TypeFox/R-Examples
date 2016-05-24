require(adimpro)
if(Sys.getenv("ImageMagick")==""){
cat("Please install ImageMagick\n")
} else {
#
#      Color space transformations
#
readline("The effect of different RGB color spaces:")

# There are different RGB color spaces, output on the screen is usually in sRGB
# while Adobe, wGamut and kodak allow to store a wider range of colors 
X11(width=10.4,height=4)
par(mfrow=c(1,4),mar=c(1,1,1,.1),mgp=c(2,1,0))
img <- read.raw(system.file("img/rawimage.png",package="adimpro"))
show.image(img,cspace="sRGB",xaxt="n",yaxt="n")
title("Output rgb space: sRGB")
show.image(img,cspace="Adobe",xaxt="n",yaxt="n")
title("Output rgb space: Adobe")
show.image(img,cspace="wGamut",xaxt="n",yaxt="n")
title("Output rgb space: wGamut")
show.image(img,cspace="kodak",xaxt="n",yaxt="n")
title("Output rgb space: kodak")

readline("Different color representations:")

# The color spaces XYZ, YUV and YIQ  are usually used for color manipulations. 
# They are related to the RGB spaces by a linear transformation.
# The color space HSI codes colors by Hue, color Saturation and color Intensity. # In package adimpro images can be converted into and from these color spaces.
X11(width=10.4,height=4)
par(mfrow=c(1,4),mar=c(1,1,1,.1),mgp=c(2,1,0))
show.image(img,cspace="xyz",xaxt="n",yaxt="n")
title("Output color space: xyz")
show.image(img,cspace="hsi",xaxt="n",yaxt="n")
title("Output color space: hsi")
show.image(img,cspace="yuv",xaxt="n",yaxt="n")
title("Output color space: yuv")
show.image(img,cspace="yiq",xaxt="n",yaxt="n")
title("Output color space: yiq")

# Function plot converts images automatically to an RGB space and,
# if neccessary, applies a gamma correction for display. 
# Function plot provides histograms of colors channels in the actual color space.
readline("Plot original Adobe RGB:")

X11(width=10,height=6)
plot(img)
readline("Plot image in sRGB RGB:")

plot(img,cspace="sRGB")
readline("Plot image in wGamut RGB:")

plot(img,cspace="wGamut")
readline("Plot image in sRGB with gamma correction:")

plot(img,cspace="sRGB",gammatype="ITU")
readline("Plot image in HSI:")

plot(rgb2hsi(img))

readline("Converting colors:")

# Function adjust.image allows for quite general color transformations,
# including changing the color space, black point, white point, 
# color temperature, exposure and gamma correction
X11(width=10.4,height=4)
par(mfrow=c(1,4),mar=c(1,1,1,.1),mgp=c(2,1,0))
show.image(img,xaxt="n",yaxt="n")
title("Original (wp D65)")
show.image(adjust.image(img,whitep="D50"),xaxt="n",yaxt="n")
title("whitepoint D50")
show.image(adjust.image(img,exposure=1.8),xaxt="n",yaxt="n")
title("increased exposure")
show.image(adjust.image(img,black=-.04,temp=6000),xaxt="n",yaxt="n")
title("changed black point and temperature")
z <- readline("End of demo, remove created objects (Y/N) :")

graphics.off()
if(toupper(z)!="N"){
rm(img, z)
}
}