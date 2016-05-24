
# create a matrix for the tests - values 0, 1 and 2 in 3x3 blocks of 10x5 numbers 012 120 201 RGB GBR BRG
ColorMatrix <- rbind(matrix(rep(0:2, each = 10, times = 5), ncol = 30, byrow = TRUE), matrix(rep(c(1, 2, 0), each = 10, times = 5), ncol = 30, byrow = TRUE), matrix(rep(c(2, 0, 1), 
    each = 10, times = 5), ncol = 30, byrow = TRUE))
GpMatrixr2gnu(ColorMatrix, "3colors.dat")  # save the matrix

# 2D heatmap with data from a file, colorbox with color gradient
Gprun("#set terminal png;set output \"2Dheatmapv1.png\"\nset view map; set size square;set yrange [] reverse\nset palette defined (0 \"red\", 1 \"green\",  2 \"blue\" )\nsplot [0:30][0:15] \"3colors.dat\" matrix with image notitle", 
    TRUE)


# 2D heatmap, the colors blend when there is a color transition, colorbox with color gradient
Gprun("#set terminal png;set output \"2Dheatmapv2.png\"\nset view map; set size square;set yrange [] reverse\nset palette defined (0 \"red\", 1 \"green\",  2 \"blue\" )\nsplot [0:30][0:15] \"3colors.dat\" matrix with pm3d notitle", 
    TRUE)


# colorbox with color gradient and plot with solid colors
Gprun("#set terminal png;set output \"2Dheatmapv3.png\"\nset view map; set size square;set yrange [] reverse\nset palette defined (0 \"red\", 1 \"green\",  2 \"blue\" )\nset pm3d corners2color c1\nsplot [0:30][0:15] \"3colors.dat\" matrix w pm3d notitle", 
    TRUE)


# Defining the palette through a file, colorbox with color gradient and plot with solid colors
Gprun("set view map; set size square;set yrange [] reverse\nset palette model RGB file \"-\"\n1 0 0\n0 1 0\n0 0 1\ne\nset pm3d corners2color c1\nsplot [0:30][0:15] \"3colors.dat\" matrix w pm3d notitle", 
    TRUE)


# colorbox with solid colors, both plot and colorbox show solid colors
Gprun("#set terminal png;set output \"2Dheatmapv4.png\"\nset view map; set size square;set yrange [] reverse\nset palette defined ( -1 \"red\", 0 \"red\", 0 \"green\", 1 \"green\", 1 \"blue\", 2 \"blue\" )\nset pm3d corners2color c1\nsplot [0:30][0:15] \"3colors.dat\" matrix w pm3d notitle", 
    TRUE)


# Defining the palette through a file, both plot and colorbox show solid colors
Gprun("set view map; set size square;set yrange [] reverse\nset palette model RGB file \"-\"\n-1 1 0 0\n0 1 0 0\n0 0 1 0\n1 0 1 0\n1 0 0 1\n2 0 0 1\ne\nset pm3d corners2color c1\nsplot [0:30][0:15] \"3colors.dat\" matrix w pm3d notitle", 
    TRUE)


# saving the palette with color gradient
Gpmatrix2palette(c(1, 0, 0, 0, 1, 0, 0, 0, 1), "3colors.pal")
# saving the palette with solid colors
Gpmatrix2palette(c(1, 0, 0, 0, 1, 0, 0, 0, 1), "3colorsSolid.pal", SolidColor = TRUE)

# Hide the colorbox and fill the gaps at the edges of the plot
GpMatrixfilePad("3colors.dat", "3colorsb.dat")  #duplicates the last column and row on the matrix data file
Gprun("#set terminal png;set output \"2Dheatmapv7.png\"\nset view map; set size square;set yrange [] reverse\nunset colorbox\nset palette model RGB file \"3colors.pal\"\nset pm3d corners2color c1\nsplot [0:30][0:15] \"3colorsb.dat\" matrix w pm3d notitle", 
    TRUE)


# create a random matrix and plot it
n <- 15
z <- matrix(sample(0:n, n^2, replace = TRUE), n, n, byrow = FALSE)
z
GpMatrixr2gnu(z, "randomMtrx.dat")
GpMatrixfilePad("randomMtrx.dat", "randomMtrxb.dat")
gpcolornames <- GpShowPaletteColornames()
randndx <- sample(1:dim(gpcolornames)[1], (n + 1))
randcolor <- as.matrix(gpcolornames[randndx, c("R", "G", "B")])
Gpmatrix2palette(c(t(randcolor)), "nrandcolorsSolid.pal", SolidColor = TRUE)
Gprun("#set terminal png;set output \"2Dheatmapv8.png\"\nset view map; set size square;set yrange [] reverse\nset palette model RGB file \"nrandcolorsSolid.pal\" u ($4):($1/255):($2/255):($3/255)\nset pm3d corners2color c3\nsplot [0:" %s% 
    n %s% "][0:" %s% n %s% "] \"randomMtrxb.dat\" matrix w pm3d notitle", TRUE)


# creating a matrix with the full color RGB values of blutux.png
gpfile <- system.file("extdata/blutux.rgb", package = "Rgnuplot")
if (!file.exists("blutux.rgb")) file.copy(gpfile, getwd())
if (!file.exists("blutux.png")) GpRGB2image("blutux.rgb", "blutux.png", 128, 128)  # convert from raw RGB to PNG

PNGdata <- GpPNG2color("blutux.png")  #get the color matrix from a PNG file
str(PNGdata)
GpMatrixr2gnu(PNGdata, "blutuxRGBfullColorMatrix.dat")  #save the color matrix to a file
GpMatrix2XYdata("blutuxRGBfullColorMatrix.dat", "blutuxRGBfullColorXYdata.dat")  #create an XY file from the color matrix file

# plot an image stored in a text datafile in XY format, with full color RGB - with points
Gprun("#set terminal pngcairo;set output \"blutuxRGBfullColor1.png\"\nset view map; set size square;set yrange [] reverse\nplot [1:128][1:128] \"blutuxRGBfullColorXYdata.dat\" u 2:1:3:($3) w p pt 5 ps 1 lc rgb var notit", 
    TRUE)


# plot an image stored in a text datafile in XY format, with full color RGB - with image
Gprun("#set terminal pngcairo;set output \"blutuxRGBfullColor2.png\"\nreset; set size square;set yrange [] reverse\nunset key;unset border;unset xtics;unset ytics\nplot [1:128][1:128] \"blutuxRGBfullColorXYdata.dat\" u 2:1:($3/65536):((int($3) & 65280)/256):(int($3) & 255) w rgbimage notit", 
    TRUE)

# plot an image stored in a text datafile in matrix format, with full color RGB
Gprun("set size square;set yrange [] reverse;\nplot \"blutuxRGBfullColorMatrix.dat\" u 1:2:($3/65536):((int($3) & 65280)/256):(int($3) & 255) matrix w rgbimage notit", TRUE)


# plot an image stored in a text datafile in matrix format, with the default palette
Gprun("#set terminal pngcairo;set output \"blutuxRGBfullColor3.png\"\nset size square;set yrange [] reverse;plot \"blutuxRGBfullColorMatrix.dat\" matrix w image notit", TRUE)

gpfile <- system.file("extdata/blutuxwithpalette.png", package = "Rgnuplot")
if (!file.exists("blutuxwithpalette.png")) file.copy(gpfile, getwd())
PNGdata256 <- GpPNG2color("blutuxwithpalette.png")  #get the color matrix from an indexed PNG file
unique(c(PNGdata256))  #236
paletteRGB <- GpcreatePaletteFromMatrix(PNGdata256)  #create a palette
paletteRGB <- c(paletteRGB, rep(0, 19))
if (file.exists("blutuxwithpalette.pal")) file.remove("blutuxwithpalette.pal")
GpRGB1to3channels(paletteRGB, fileRGB3channel = "blutuxwithpalette.pal")  #save the palette to a file with separated RGB components
GpPalettePlot("blutuxwithpalette.png", "GIMP")  #show the palette from the PNG file
PNGdataIndexed <- GpcreateIndexFromMatrixAndPalette(PNGdata256, paletteRGB)  #create an indexed color matrix
GpMatrixr2gnu(PNGdataIndexed, "blutuxwithpalette.dat")  #save the indexed color matrix to a file
# plot an image stored in a text datafile in matrix format, with its own palette
max(PNGdataIndexed)
Gprun("#set terminal pngcairo;set output \"blutuxwithpalette1.png\"\nset view map; set size square;set yrange [] reverse\nset cbrange [0:254]\nset palette model RGB file \"blutuxwithpalette.pal\" u ($1/255):($2/255):($3/255)\nsplot \"blutuxwithpalette.dat\" matrix w image notit", 
    TRUE)


if (!file.exists("blutux.png")) GpRGB2image(system.file("extdata/blutux.rgb", package = "Rgnuplot"), "blutux.png", 128, 128)
testmap <- "blutux"  #tuxgnu file to be mapped
# GpimageResize(testmap %s% '0.png',testmap %s% '256.png',256,128)
PNGdata <- GpPNG2color(testmap %s% ".png")  #get the color matrix
unique(c(PNGdata))  # 418 different colors
paletteRGBx <- GpcreatePaletteFromMatrix(PNGdata)  #create a palette
GpRGB1to3channels(paletteRGBx, fileRGB3channel = "blutux.pal")  #save the palette to a file with separated RGB components
PNGdataIndexedx <- GpcreateIndexFromMatrixAndPalette(PNGdata, paletteRGBx)  #create an indexed color matrix
# str(PNGdataIndexed)
GpMatrixr2gnu(PNGdataIndexedx, "blutux.dat")  #save the indexed color matrix to a file
Gprun("#set terminal pngcairo;set output \"blutuxwithpalette2.png\"\nset view map; set size square;set yrange [] reverse\nset cbrange [0:" %s% (length(unique(c(PNGdata))) - 1) %s% "]\nset palette model RGB file \"blutux.pal\" u ($1/255):($2/255):($3/255)\nsplot \"blutux.dat\" matrix w image notit", 
    TRUE)


# using a GIMP palette on an indexed PNG, color effects are possible with cbrange
Gprun("set view map;set size ratio -1;set lmargin 0;set rmargin 0;set tmargin 0;set bmargin 0;unset key;unset tics;unset border;set yrange [] reverse\nset palette model RGB file \"blutuxwithpalettePNG.gpl\" every ::5 u ($1/255):($2/255):($3/255)\nset cbrange [0:255]\nsplot \"blutuxwithpalette.dat\" matrix w image notit\n", 
    TRUE)

# splot 'blutuxwithpalette.png' binary filetype=png w rgbimage notitle


# using a palette on an indexed PNG with pm3d, notice the complexity of the palette on the colorbox
Gprun("#set terminal pngcairo;set output \"blutuxwithpalette3.png\"\nset view map;\nset size ratio -1;set lmargin 0;set rmargin 0;set tmargin 0;set bmargin 0;unset key;unset tics;unset border;set yrange [] reverse\nset palette model RGB file \"blutuxwithpalette.pal\" u ($1/255):($2/255):($3/255)\nset cbrange[0:254]\nset pm3d corners2color c1\nsplot \"blutuxwithpalette.dat\" matrix with pm3d notitle", 
    TRUE)


# plot Tux with palette into a hemisphere
Gprun("#set terminal pngcairo;set output \"blutuxwithpalette4.png\"\nset angles degrees\nunset colorbox\nset parametric\nset isosamples 11,11\nset mapping spherical\nset yzeroaxis linetype 0 linewidth 1.000 \nset urange [ -90.0000 : 90.0000 ]\nset vrange [ 0.00000 : 360.000 ]\nset hidden3d back\nset view equal xyz\nset view 86,86,1,1\nset palette model RGB file \"blutuxwithpalette.pal\" u ($1/255):($2/255):($3/255)\nset cbrange[0:254]\nset pm3d corners2color c1\nsplot cos(u)*cos(v),cos(u)*sin(v),sin(u) w l lc rgb \"grey\" notit, \"blutuxwithpalette.dat\" matrix u (($1-64)/128*180):(-($2-64)/128*180):(1):3 w pm3d notit", 
    TRUE)

# plot Tux with palette into a hemisphere XY file
GpMatrix2XYdata("blutuxwithpalette.dat", "blutuxwithpaletteXY.dat")
Gprun("set angles degrees\nunset colorbox\nunset tics\nunset border\nset parametric\nset isosamples 11,11\nset mapping spherical\nset yzeroaxis linetype 0 linewidth 1.000 \nset urange [ -90.0000 : 90.0000 ]\nset vrange [ 0.00000 : 360.000 ]\nset hidden3d back\nset view equal xyz\nset view 80,80,2,2\nset palette model RGB file \"blutuxwithpalette.pal\" u ($1/255):($2/255):($3/255)\nset cbrange[0:254]\nset pm3d corners2color c1\nsplot cos(u)*cos(v),cos(u)*sin(v),sin(u) w l lc rgb \"grey\" notit, \"blutuxwithpaletteXY.dat\" u (($2-64)/128*180):(-($1-64)/128*180):(1):3 w pm3d notit", 
    TRUE)



# plotting an RGB image as RGB or RGBA returns the same output because the alpha channel is considered zero for all pixels, all colors are opaque.

# download an image of the GNU mascot from Wikipedia
gnuimg <- "500px-Heckert_GNU_white.png"
download.file("http://upload.wikimedia.org/wikipedia/en/thumb/2/22/Heckert_GNU_white.svg/500px-Heckert_GNU_white.svg.png", gnuimg)
# the image has an alpha channel, if ignored then the background on the transparent area becomes black
GpImagePlot(gnuimg)
# if the alpha channel is removed then the background on the transparent area becomes white
GpImagePlot(gnuimg, TRUE)
# plot the image with a red background
GpImagePlot(gnuimg, TRUE, "red")
# create a new image with green background
GpImage2PNG(gnuimg, "gnugreen.png", FALSE, TRUE, "green")
# resize the image to 128x128
GpimageResize("gnugreen.png", "gnu.png", 128, 128)
# tile 2 images
GpImageTile("tuxgnu.png", matrix(c("blutux.png", "gnu.png"), 2, 1), c(128, 128), c(128))
GpImagePlot("tuxgnu.png")

testmap <- "tuxgnu"  # file to be mapped
# GpimageResize(testmap %s% '0.png',testmap %s% '256.png',256,128)
PNGdata <- GpPNG2color(testmap %s% ".png")  #get the color matrix
paletteRGBx <- GpcreatePaletteFromMatrix(PNGdata)  #create a palette
NpaletteColors <- length(paletteRGBx) - 1  # number of palette colors starting from zero
GpRGB1to3channels(paletteRGBx, fileRGB3channel = testmap %s% ".pal")  #save the palette to a file with separated RGB components
PNGdataIndexedx <- GpcreateIndexFromMatrixAndPalette(PNGdata, paletteRGBx)  #create an indexed color matrix
# str(PNGdataIndexed)
GpMatrixr2gnu(PNGdataIndexedx, testmap %s% ".dat")  #save the indexed color matrix to a file

# plot 6 views of the sphere
Gprun("set angles degrees\nunset colorbox\nunset border\nunset xtics\nunset ytics\nunset ztics\nset parametric\nset isosamples 11,11\nset mapping spherical\nset yzeroaxis linetype 0 linewidth 1.000 \n#set hidden3d back\nset view equal xyz\nset palette model RGB file \"tuxgnu.pal\" u ($1/255):($2/255):($3/255)\nset cbrange[0:" %s% 
    NpaletteColors %s% "]\nset pm3d corners2color c1\nset multiplot layout 2,3\nset view 360,0,2,2\nsplot \"tuxgnu.dat\" matrix u (($1-64)/128*180):(-($2-64)/128*180):(1):3 w pm3d notit\nset view 180,0,2,2\nsplot \"\" matrix u (($1-64)/128*180):(-($2-64)/128*180):(1):3 w pm3d notit\nset view 91,90,2,2\nsplot \"\" matrix u (($1-64)/128*180):(-($2-64)/128*180):(1):3 w pm3d notit\nset view 89,270,2,2\nsplot \"\" matrix u (($1-64)/128*180):(-($2-64)/128*180):(1):3 w pm3d notit\nset view 93,0,2,2 #87,0,2,2\nsplot \"\" matrix u (($1-64)/128*180):(-($2-64)/128*180):(1):3 w pm3d notit\nset view 87,180,2,2\nsplot \"\" matrix u (($1-64)/128*180):(-($2-64)/128*180):(1):3 w pm3d notit\nunset multiplot", 
    TRUE)
 
