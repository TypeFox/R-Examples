
# download an image of a windrose from Wikipedia
download.file("http://upload.wikimedia.org/wikipedia/commons/thumb/1/1a/Brosen_windrose.svg/500px-Brosen_windrose.svg.png", "windrose.png")
# the image seems to have a black background, however the background is transparent because this image is RGBA (RGB+alpha channel)
GpImagePlot("windrose.png")
# the 'alpha' parameter will show the image correctly
GpImagePlot("windrose.png", TRUE)
# by converting from RGBA to RGB, the transparent background is replaced by a white one
GpImage2PNG("windrose.png", "windrose2.png", alpha = TRUE)
# now the background is white
GpImagePlot("windrose2.png")

if (!file.exists("blutux.png")) {
    file.copy(system.file("extdata/blutux.rgb", package = "Rgnuplot"), getwd())
    GpRGB2image("blutux.rgb", "blutux.png", 128, 128)
}
GpImageRgbfiltercolorSepia("blutux.png", "blutuxsepia.png")
GpImageRgbgreyscaleBT709("blutux.png", "blutuxgreyscaleBT709.png")
GpImageRgbfalsecolor("blutux.png", "blutuxfalsecolor.png")
# make an image lighter by increasing all of the RGB channel values by 50%, limit the values to 255
GpimageRGBchange("blutux.png", "blutuxlighter.png", "PNG", NULL, "(1.5*r>255)?255:(1.5*r)", "(1.5*g>255)?255:(1.5*g)", "(1.5*b>255)?255:(1.5*b)")
# create an image from tiling other images
GpImageTile("tux4x4.png", matrix(c("blutuxsepia.png", "blutuxgreyscaleBT709.png", "blutuxfalsecolor.png", "blutuxlighter.png"), 2, 2), c(128, 128), c(128, 128))
GpImagePlot("tux4x4.png")  # GpImagePlot('tux4x4.png',filetype='PNG') 
