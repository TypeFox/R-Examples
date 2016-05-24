library(maps)
library(mapdata)
# if it doesn't exist, create a file with the volcano data
if (!file.exists("volcano.txt")) write.table(volcano, file = "volcano.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 3D wireframe
Gprun("#set terminal postscript eps color;set output \"Maungawhau1.eps\"\nset hidden3d\nsplot \"volcano.txt\" matrix w l", TRUE)

# 3D color surface
Gprun("#set terminal postscript eps color;set output \"Maungawhau2.eps\"\nset pm3d\nsplot \"volcano.txt\" matrix w pm3d", TRUE)

# 2D heatmap
Gprun("#set terminal postscript eps color;set output \"Maungawhau3.eps\"\nset pm3d map\nsplot \"volcano.txt\" matrix w pm3d", TRUE)

# 2D countour plot
Gprun("#set terminal png;set output \"Maungawhau4.png\"\n#unset clabel;\nunset surface\nset view map\nset contour base\nsplot \"volcano.txt\" u 2:1:3 matrix w l notitle", TRUE)

# 3D wireframe + 2D heatmap
Gprun("#set terminal png;set output \"Maungawhau5.png\"\nset pm3d at b\nsplot \"volcano.txt\" matrix w l", TRUE)

# 3D color surface + 2D contour
Gprun("#set terminal png;set output \"Maungawhau6.png\"\nset contour base\nset hidden3d\nsplot \"volcano.txt\" matrix w pm3d nocontour, \"\" matrix  w lines nosurface", TRUE)

# 3D surface + contour
Gprun("#set terminal postscript eps color;set output \"Maungawhau7.eps\"\nset contour surface\nset hidden3d\nsplot \"volcano.txt\" matrix w pm3d notitle", TRUE)

# 3D surface + black contour
tmpfile1 <- tempfile()
tmpfile2 <- tempfile()
Gprun("#set terminal postscript eps color;set output \"Maungawhau8.eps\"\nset contour base; set cntrparam level 5\nunset surface\nunset clabel\nset table \"" %s% tmpfile1 %s% "\"\nsplot \"volcano.txt\" matrix w l lt -1\nunset table\nreset\n!awk \"NF<2{printf\\\"\\n\\\"}{print}\" < " %s% 
    tmpfile1 %s% " > " %s% tmpfile2 %s% "\nsplot \"volcano.txt\" matrix w pm3d nocontour notit,\"" %s% tmpfile2 %s% "\" w l lc rgb \"black\" nosurface notit", TRUE)

# satellite image of Maunga Whau
library(RgoogleMaps)
# save the Google Maps tile
m <- GetMap(center = c(-36.875673, 174.765036), zoom = 16, destfile = "MaungawhauZOOM16.png", maptype = "satellite")
GpImagePlot("MaungawhauZOOM16.png")

# use GIMP to crop and rescale crop 303x434
GpImageCrop("MaungawhauZOOM16.png", "MaungawhauZOOM16crop.png", 150, 92, 453, 526)

# rotate 180 degrees
GpImageRotate("MaungawhauZOOM16crop.png", "MaungawhauZOOM16cropRot.png", 180)

# resize to 61x87
GpimageResize("MaungawhauZOOM16cropRot.png", "Maungawhau61x87.png", 61, 87)

GpImagePlot("Maungawhau61x87.png")

# create a gnuplot data file with color and DEM
filePNG <- "Maungawhau61x87.png"  #system.file('extdata/Maungawhau2.png',package='Rgnuplot')
fileDEMtab <- "volcano.txt"
file3Ddat <- "Maungawhau2.dat"  # color and DEM data file
GpPNG4DEM(filePNG, fileDEMtab, file3Ddat)  # create the data file

# 3D DEM with points
Gprun("#set terminal png;set output \"Maungawhau10.png\"\nset view equal_axes xy\nsplot \"Maungawhau2.dat\" u 1:2:3:4 w p pt 5 ps 1 lc rgb var notitle", TRUE)

# 3D DEM with palette
Gprun("#set terminal png;set output \"Maungawhau11.png\"\nset palette defined (0 \"black\",1 \"dark-green\", 2 \"dark-olivegreen\",3 \"grey50\", 4 \"white\")\nset view equal_axes xy\nunset colorbox\nset hidden3d\nsplot \"Maungawhau2.dat\" u 1:2:3:4 w pm3d notitle", 
    TRUE)

# create a palette for Maungawhau61x87.png
testmap <- "Maungawhau61x87"
PNGdata <- GpPNG2color(testmap %s% ".png")  # get the color matrix from the PNG file
paletteRGB <- GpcreatePaletteFromMatrix(PNGdata)  # create a palette
GpRGB1to3channels(paletteRGB, fileRGB3channel = testmap %s% ".pal")  # save the palette to a file with separated RGB components
# create an indexed color matrix for Maungawhau61x86.png
PNGdataIndexed <- GpcreateIndexFromMatrixAndPalette(PNGdata, paletteRGB)  # create an indexed color matrix
p <- matrix(c(t(PNGdataIndexed)), ncol = 1)
GpMatrix2XYdata("volcano.txt", testmap %s% "XYndx.dat", p, surfacegrid = FALSE, TRUE)
NpaletteColors <- length(paletteRGB) - 1  # number of palette colors starting from zero

Gprun("#set terminal png;set output \"Maungawhau12.png\"\nset palette model RGB file \"" %s% testmap %s% ".pal\" u ($1/255):($2/255):($3/255)\nset cbrange[0:" %s% NpaletteColors %s% 
    "]\nset view equal xyz\nset view 38,26\nset pm3d corners2color c1\nunset colorbox\nset hidden3d\nsplot \"" %s% testmap %s% "XYndx.dat\" u 1:2:3:4 w pm3d notitle", TRUE)

## download a Google Maps tile and use its coordinates to clip a DEM (in geoTIFF format) New Zealand National Digital Elevation Model (North Island) from Landcare Research, Sourced
## from LINZ Geodetic Database. Crown Copyright Reserved.  from http://lris.scinfo.org.nz/layer/131-nzdem-north-island-25-metre/ save data from coords: 174.7581, 174.7718, -36.88127,
## -36.87013 the file saved as 1-0001-0001.tif

# save the Google Maps tile (if it was not already done before)
m1 <- GetMap(center = c(-36.875673, 174.765036), zoom = 16, destfile = "MaungawhauZOOM16.png", maptype = "satellite")

# Read the geoTiff and save it to a data file with the coordinates only and without the coordinate reference system (CRS) arguments.
f <- "1-0001-0001.tif"
if (file.exists(f)) {
    require(rgdal)
    x <- readGDAL(f)
    x2 <- raster(x)
    bbx = m1$BBOX
    r <- crop(x2, extent(bbx$ll[, 2], bbx$ur[, 2], bbx$ll[, 1], bbx$ur[, 1]))
    r2 <- matrix(r@data@values, dim(r)[1], dim(r)[2], byrow = TRUE)
    # r2<-apply(r2,2,rev) #flip the matrix
    GpSplot(r2)
    str(r)  #44 rows x 54 columns
    GpMatrixr2gnu(r2, "Maungaw.dat")
    Gprun("splot \"Maungaw.dat\" matrix w l", TRUE)
}
# Instead of downloading the geoTiff and running code , the file 'Maungaw.dat' can be copied from the extdata directory of Rgnuplot.
if (!(file.exists(f))) file.copy(system.file("extdata/Maungaw.dat", package = "Rgnuplot"), getwd())

# resample to five times its original size.
GpResampleDEM("Maungaw.dat", "MaungawresXY.dat", c(0, 44 * 5 - 1), c(0, 54 * 5 - 1), interpolationMethod = (44 * 5) %s% "," %s% (54 * 5) %s% " gauss 1")
Gprun("splot \"MaungawresXY.dat\" u 1:2:3 w l", TRUE)

# Rotate and resize the satellite image.
GpImageRotate("MaungawhauZOOM16.png", "MaungawhauZOOM16rotated.png", 270)
GpimageResize("MaungawhauZOOM16rotated.png", "Maungawhau220x270.png", 44 * 5, 54 * 5)
GpImagePlot("Maungawhau220x270.png")

pngmap <- "Maungawhau220x270"
PNGdata <- GpPNG2color(pngmap %s% ".png")  # get the color matrix from the PNG file
paletteRGB <- GpcreatePaletteFromMatrix(PNGdata)  # create a palette
GpRGB1to3channels(paletteRGB, fileRGB3channel = pngmap %s% ".pal")  # save the palette to a file with separated RGB components
PNGdataIndexed <- GpcreateIndexFromMatrixAndPalette(PNGdata, paletteRGB)  # create an indexed color matrix
NpaletteColors <- length(paletteRGB) - 1  # number of palette colors starting from zero

# create a gnuplot data file with color and DEM
GpXydata2matrix("MaungawresXY.dat", "Maungaw640.dat")
Gprun("set pm3d;splot \"Maungaw640.dat\" matrix w pm3d", TRUE)
p <- matrix(c(t(PNGdataIndexed)), ncol = 1)
GpMatrix2XYdata("Maungaw640.dat", "Maungawhau3.dat", p, surfacegrid = FALSE, TRUE)

Gprun("#set terminal pngcairo;set output \"Maungawhau13.png\"\nset palette model RGB file \"" %s% pngmap %s% ".pal\" u ($1/255):($2/255):($3/255)\nset cbrange[0:" %s% NpaletteColors %s% 
    "]\n#set view equal xyz\nunset tics;unset border\nset size ratio 1,1,5\nset view 30,130\nset pm3d corners2color c1\nunset colorbox\nset hidden3d\nsplot \"Maungawhau3.dat\" u 1:2:3:4 w pm3d notitle", 
    TRUE)
 
