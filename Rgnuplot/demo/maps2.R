library(maps)
library(mapdata)
# saving a map from package 'maps' with Mercator projection
map("world", projection = "mercator")
if (!file.exists("worldRmapMercator.dat")) {
    GpMapsr2gnu(map("world", projection = "mercator", plot = FALSE), "worldRmapMercator.dat")
}
Gprun("unset key; unset tics\nplot \"worldRmapMercator.dat\" w l lc rgb \"black\"", TRUE)

if (!file.exists("ne_110m_coastline.dat")) {
    # converting a world map from SHP to a data file readable by gnuplot download from http://www.naturalearthdata.com/downloads/ ne_110m_coastline.zip only 84 kilobytes
    if (!file.exists("ne_110m_coastline.zip")) 
        download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/physical/ne_110m_coastline.zip", "ne_110m_coastline.zip", method = "wget")
    if (file.info("ne_110m_coastline.zip")$size == 0) 
        stop("Please download http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/physical/ne_110m_coastline.zip")
    unzip("ne_110m_coastline.zip", exdir = "ne_110m_coastline")
    world.shp <- readShapeLines("ne_110m_coastline/ne_110m_coastline.shp", proj4string = CRS("+proj=longlat"))
    GpSHP2gnu("ne_110m_coastline", "ne_110m_coastline", "ne_110m_coastline.dat")  # it's already WGS84
    plot(world.shp)
}
Gprun("#set term png;set output \"naturalearthdata_1.png\"\nplot \"ne_110m_coastline.dat\" w l lc rgb \"black\"", TRUE)

if (!file.exists("worldpar.dat")) GpMapMerpar("worldpar.dat", "worldmer.dat", 10, 10)
Gprun("#set term png;set output \"naturalearthdata_2.png\"\nunset key; unset tics;unset border\nplot \"ne_110m_coastline.dat\" w l lc rgb \"black\", \"worldpar.dat\"  w l lc rgb \"grey\", \"worldmer.dat\"  w l lc rgb \"grey\"", 
    TRUE)

# convert the shapefile with Tissot's indicatrices to a data file format using the shapefile from Matthew T. Perry 2005, Tissot Indicatrix - Examining the distortion of map
# projections http://blog.perrygeo.net/2005/12/11/tissot-indicatrix-examining-the-distortion-of-2d-maps/
tissfiles <- dir(system.file(package = "Rgnuplot") %s% "/extdata", pattern = "tissot")
if (!all(file.exists(tissfiles))) file.copy(system.file(paste("extdata/", tissfiles, sep = ""), package = "Rgnuplot"), getwd())
GpSHP2gnu(".", "tissot", "tissot.dat")

if (!file.exists("NOAACoastline.dat")) {
    download.file("http://www.ngdc.noaa.gov/mgg/coast/tmp/13777.dat", "NOAA Coastline Data.dat", method = "wget")
    if (!file.exists("NOAA Coastline Data.dat")) 
        stop("Please download WCL(World Coast Line) from http://www.ngdc.noaa.gov/mgg_coastline/ and save it as \"NOAA Coastline Data.dat\"")
    Gprun("!awk '{gsub(/>/,\"\")}1' \"NOAA Coastline Data.dat\" >NOAACoastline.dat", FALSE)
    
}

Gprun("unset key;plot \"NOAACoastline.dat\" w l, \"worldpar.dat\" w l,\"worldmer.dat\" w l", TRUE)

# Mercator's Projection
Gprun("#set term png;set output \"NOAACoastline_Mercator.png\"\nload \"projections.gnu\"\np=MercatorInit(0)\nset size ratio -1\nset title p\nplot \"worldmer.dat\" using (MercatorYC($2,$1)):(MercatorXC($2,$1)) notit w l ls 3, \\\n \"worldpar.dat\" using (MercatorYC($2,$1)):(MercatorXC($2,$1)) notit w l ls 2, \\\n \"NOAACoastline.dat\" using (MercatorYC($2,$1)):(MercatorXC($2,$1)) notit w l ls 1", 
    TRUE)

# Mercator's Projection with Tissot's indicatrices
Gprun("#set term png;set output \"NOAACoastline_Mercator_Tissot.png\"\nload \"projections.gnu\"\np=MercatorInit(0)\nset size ratio -1\nset title p\nplot \"tissot.dat\" using (MercatorYC($2,$1)):(MercatorXC($2,$1)) w l notit, \\\n\"worldmer.dat\" using (MercatorYC($2,$1)):(MercatorXC($2,$1)) notit w l ls 3, \\\n \"worldpar.dat\" using (MercatorYC($2,$1)):(MercatorXC($2,$1)) notit w l ls 2, \\\n \"NOAACoastline.dat\" using (MercatorYC($2,$1)):(MercatorXC($2,$1)) notit w l ls 1", 
    TRUE)


# download the Cloudless Earth (Day) image file from NASA
if (!file.exists("earth_day.jpg")) download.file("http://nssdc.gsfc.nasa.gov/planetary/image/earth_day.jpg", "earth_day.jpg")
testmap <- "earth_day"
if (!file.exists(testmap %s% ".png")) {
    print("Please wait, this might take some time.")
    # GpImage2PNG(testmap %s% '.jpg',testmap %s% '.png') # convert from JPEG to PNG
    GpimageResize(testmap %s% ".jpg", testmap %s% ".png", 500, 250)  #reduce the size to 1/4 (1000,500  1/2 250,125 1/8)
}
PNGdata2 <- GpPNG2color(testmap %s% ".png")  #get the color matrix from the PNG file
Mheight <- dim(PNGdata2)[1]
Mwidth <- dim(PNGdata2)[2]
paletteRGB <- GpcreatePaletteFromMatrix(PNGdata2)  #create a palette
NpaletteColors <- length(paletteRGB) - 1  # number of palette colors starting from zero
if (!file.exists(testmap %s% ".pal")) {
    GpRGB1to3channels(paletteRGB, fileRGB3channel = testmap %s% ".pal")  # save the palette to a file with separated RGB components
}
PNGdataIndexed <- GpcreateIndexFromMatrixAndPalette(PNGdata2, paletteRGB)  # create an indexed color matrix
if (!file.exists(testmap %s% "matrixndx.dat")) {
    print("Please wait, this might take some time.")
    GpMatrixr2gnu(PNGdataIndexed, testmap %s% "matrixndx.dat")  # save the indexed color matrix to a file
}
NpaletteColors <- length(paletteRGB) - 1  # number of palette colors starting from zero
# using splot to plot a 2D map from a data file in XY format
if (!file.exists(testmap %s% "XYndx.dat")) {
    GpMatrix2XYdata(testmap %s% "matrixndx.dat", testmap %s% "XYndx.dat")  # create an XY file from the matrix file 
}
if (!file.exists("earth_dayXYcoords.dat")) {
    # convert the coordinates and save to a new file
    GpXYcoordsConvertFun("earth_dayXYndx.dat", "earth_dayXYcoords.dat", function(y) -(y/Mheight * 180 - 90), function(x) (x/Mwidth * 360 - 180), TRUE)
}

# create meridian and parallel lines for splot
if (!file.exists("worldparS15.dat")) GpMapMerpar("worldparS15.dat", "worldmerS15.dat", 15, 15, TRUE)
if (!file.exists("worldparS20.dat")) GpMapMerpar("worldparS20.dat", "worldmerS20.dat", 20, 20, TRUE)
# Mercator's Projection
Gprun("#set term pngcairo;set output \"earth_dayMercator.png\"\nload \"projections.gnu\"\np=MercatorInit(0)\nset title p\nset pm3d map\nunset key;unset tics;unset border;unset colorbox\nset palette model RGB file \"" %s% 
    testmap %s% ".pal\" u ($1/255):($2/255):($3/255)\nset cbrange[0:" %s% NpaletteColors %s% "]\nset pm3d corners2color c1\nsplot \"earth_dayXYcoords.dat\" u (MercatorYC($2,$1)):(MercatorXC($2,$1)):3 w pm3d notit, \\\n\"worldmerS15.dat\" using (MercatorYC($2,$1)):(MercatorXC($2,$1)):(1) notit w l ls 3, \\\n\"worldparS15.dat\" using (MercatorYC($2,$1)):(MercatorXC($2,$1)):(1) notit w l ls 2", 
    TRUE)

 
