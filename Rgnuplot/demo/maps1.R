library(maps)
library(mapdata)
if (!file.exists("worldRmap.dat")) {
    # convert the map of the world from the R package maps into a format readable by gnuplot
    GpMapsr2gnu(map("world", plot = FALSE), "worldRmap.dat")
}
# plot the map of the world - cartesian coordinate system by default - Equirectangular Projection (plate carree projection)
Gprun("#set terminal png;set output \"worldRmap.png\"\nunset key\nunset xtics\nunset ytics\nset xrange[-179:179]#longitude\nset yrange[-89:89]#latitude\nplot \"worldRmap.dat\" w l lc rgb \"black\"", 
    TRUE)

if (!file.exists("NewZealand.dat")) {
    # convert the map of New Zealand from the R package mapdata into a format readable by gnuplot
    GpMapsr2gnu(map("world2Hires", c("New Zealand:Stewart Island", "New Zealand:South Island", "New Zealand:North Island", "New Zealand:Great Barrier Island", "New Zealand:Auckland Island"), 
        plot = FALSE), "NewZealand.dat")
}
# plot the map of the world in black and New Zealand in red
Gprun("#set terminal png;set output \"worldNewZealand.png\"\nunset key\nunset xtics\nunset ytics\nset xrange[-179:179]#longitude\nset yrange[-89:89]#latitude\nplot \"worldRmap.dat\" w l lc rgb \"black\", \"NewZealand.dat\" w filledcurve lc rgb \"red\"", 
    TRUE)

# find the bounding box around New Zealand
NZbox <- GpBoxXY("NewZealand.dat")
# plot the map of New Zealand
Gprun("#set terminal png;set output \"NewZealandboxXY.png\"\nunset key\nunset xtics\nunset ytics\nset size ratio 1.8\nset xrange[" %s% NZbox[1, 1] %s% ":" %s% NZbox[2, 1] %s% "]#longitude\nset yrange[" %s% 
    NZbox[1, 2] %s% ":" %s% NZbox[2, 2] %s% "]#latitude\nplot \"NewZealand.dat\" w l lc rgb \"red\"", TRUE)

# plot the map of the world in black and New Zealand in red - Polar Stereographic Projection
Gprun("#set terminal png;set output \"worldNewZealandPolarStereographic.png\"\nset polar\nset angles degrees\nunset border # hide border\nunset xtics # hide X-numbers\nunset ytics # hide Y-numbers\nset hidden3d\nset grid polar 45\nset size ratio -1\nset rrange [-89:0]\nset trange [-180:180]\nset xrange[89:-89]#longitude\nset yrange[-89:89]#latitude\nset multiplot layout 1,2\nset title \"Southern Hemisphere\"\nplot (0) w l lc rgb \"black\" notit, \"worldRmap.dat\" w l lc rgb \"blue\" notit, \"NewZealand.dat\" w l lc rgb \"red\" notit\nset rtics (0,20,40,60,80)\nset rtics axis mirror \nset rrange [89:179]\nset trange [-180:180]\nset xrange[89:-89]\nset yrange[89:-89]\nset title \"Northern Hemisphere\"\nplot (0) w l lc rgb \"black\" notit, \"worldRmap.dat\" u 1:(($2>0)?$2:0/0)  w l lc rgb \"blue\" notit\nunset multiplot", 
    TRUE)

# plot the world and New Zealand in 3-dimensions - spherical - Orthographic Projection
Gprun("#set terminal png;set output \"worldNewZealandOrthographic.png\"\nset angles degrees\nunset key\nset parametric\nset samples 32 \nset isosamples 13, 25 #12 parallels 24 meridians\nset mapping spherical\nset yzeroaxis linetype 0 linewidth 1.000 \nset urange [ -90.0000 : 90.0000 ]\nset vrange [ -180.000 : 180.000 ]\nset hidden3d front\nset view equal xyz\nset view 110,250\nsplot cos(u)*cos(v),cos(u)*sin(v),sin(u) w l lc rgb \"grey\" notit, \"NewZealand.dat\" u 1:2:(1) w l lc rgb \"red\" notit,\"worldRmap.dat\" u 1:2:(1) w l lc rgb \"black\" notit", 
    TRUE)

# plot the world and New Zealand in 3-dimensions - cylindrical
Gprun("set angles degrees\nunset key\nset parametric\nset samples 32,32 \nset isosamples 13, 13# 12 parallels 12 meridians\nset mapping cylindrical\nset hidden3d back offset 0 trianglepattern 3 undefined 1 altdiagonal bentover\nset yzeroaxis linetype 0 linewidth 1.000\nset ticslevel 0\nset urange [ -180.000 : 180.000 ] noreverse nowriteback\nset vrange [ -90.0000 : 90.0000 ] noreverse nowriteback\nset view 64,250\nCyldiameter=1\nsplot cos(u)*Cyldiameter,sin(u)*Cyldiameter,v w l lc rgb \"grey\" notit, \"NewZealand.dat\" u 1:2:(Cyldiameter) w l lc rgb \"red\" notit,\"worldRmap.dat\" u 1:2:(Cyldiameter) w l lc rgb \"black\" notit", 
    TRUE)

# convert a shapefile (SHP) to a format readable by gnuplot the shapefile is the NZGD2000 Meridional Circuits, from
# http://data.linz.govt.nz/layer/817-nz-meridional-circuit-boundaries-nzgd2000/
og <- "/lds-nz-meridional-circuit-boundaries-nzgd2000-SHP"
if (!file.exists(og)) print("Please download the file \"lds-nz-meridional-circuit-boundaries-nzgd2000-SHP\"") else {
    GpSHP2gnu(og, "nz-meridional-circuit-bou", "NZ.dat")  # NZGD2000 to WGS84
    # plot the map of New Zealand with meridional circuit boundaries
    Gprun("#set terminal png;set output \"worldNewZealandMeridionalCircuit.png\"\nunset key\nunset xtics\nunset ytics\nset size ratio 1.8\nset xrange[" %s% NZbox[1, 1] %s% ":" %s% NZbox[2, 
        1] %s% "]#longitude\nset yrange[" %s% NZbox[1, 2] %s% ":" %s% NZbox[2, 2] %s% "]#latitude\nplot \"NewZealand.dat\" w l lc rgb \"red\", \"NZ.dat\" w l", TRUE)
}

# plot the map of the world - cartesian coordinate system by default, Equirectangular Projection Cloudless Earth (Day) http://nssdc.gsfc.nasa.gov/planetary/image/earth_day.jpg
# Cloudless Earth (Night) http://nssdc.gsfc.nasa.gov/planetary/image/earth_night.jpg download the Cloudless Earth (Day) image file from NASA
if (!file.exists("earth_day.jpg")) download.file("http://nssdc.gsfc.nasa.gov/planetary/image/earth_day.jpg", "earth_day.jpg")
Gprun("#set terminal pngcairo;set output \"gp_earth_day.png\"\nunset key\nunset xtics\nunset ytics\nset xrange[1:2000]\nset yrange[1:1000]\nplot \"earth_day.jpg\" binary filetype=jpg w rgbimage notit,\"worldRmap.dat\" u (($1+360/2)/360*2000):(($2+180/2)/180*1000) w l lc rgb \"red\" notit", 
    TRUE)

# create files in matrix and XY format, with the indexed colors from a PNG file, create also a palette file
testmap <- "earth_day"  # file to be mapped
if (!file.exists(testmap %s% ".png")) {
    print("Please wait, this might take some time.")
    # GpImage2PNG(testmap %s% '.jpg',testmap %s% '.png') # convert from JPEG to PNG
    GpimageResize(testmap %s% ".jpg", testmap %s% ".png", 500, 250)  #reduce the size to 1/4 (1000,500  1/2 250,125 1/8)
}
PNGdata2 <- GpPNG2color(testmap %s% ".png")  # get the color matrix from the PNG file
Mheight <- dim(PNGdata2)[1]  # map height
Mwidth <- dim(PNGdata2)[2]  # map width
paletteRGB <- GpcreatePaletteFromMatrix(PNGdata2)  # create a palette
if (!file.exists(testmap %s% ".pal")) {
    GpRGB1to3channels(paletteRGB, fileRGB3channel = testmap %s% ".pal")  # save the palette to a file with separated RGB components
}
PNGdataIndexed <- GpcreateIndexFromMatrixAndPalette(PNGdata2, paletteRGB)  # create an indexed color matrix
if (!file.exists(testmap %s% "matrixndx.dat")) {
    print("Please wait, this might take some time.")
    GpMatrixr2gnu(PNGdataIndexed, testmap %s% "matrixndx.dat")  # save the indexed color matrix to a file
}
NpaletteColors <- length(paletteRGB) - 1  # number of palette colors starting from zero

# using splot to plot a 2D map from a data file in matrix format
Gprun("set pm3d map\nset size ratio -1;set lmargin 0;set rmargin 0;set tmargin 0;set bmargin 0;unset key;unset tics;unset border;set yrange [] reverse;unset colorbox\nset palette model RGB file \"" %s% 
    testmap %s% ".pal\" u ($1/255):($2/255):($3/255)\nset cbrange[0:" %s% NpaletteColors %s% "]\nset pm3d corners2color c1\nsplot \"earth_daymatrixndx.dat\" matrix u 1:2:3 w pm3d notit", 
    TRUE)

# using splot to plot a 2D map from a data file in XY format
if (!file.exists(testmap %s% "XYndx.dat")) {
    GpMatrix2XYdata(testmap %s% "matrixndx.dat", testmap %s% "XYndx.dat")  # create an XY file from the matrix file 
}
Gprun("set pm3d map\nset size ratio -1;set lmargin 0;set rmargin 0;set tmargin 0;set bmargin 0;unset key;unset tics;unset border;unset colorbox\nset palette model RGB file \"" %s% testmap %s% 
    ".pal\" u ($1/255):($2/255):($3/255)\nset cbrange[0:" %s% NpaletteColors %s% "]\nset pm3d corners2color c1\nsplot \"earth_dayXYndx.dat\" u (($2)/1000*360-180):(-($1)/500*180-90):3 w pm3d notit", 
    TRUE)

if (!file.exists("earth_dayXYcoords.dat")) {
    # convert the coordinates and save to a new file
    GpXYcoordsConvertFun("earth_dayXYndx.dat", "earth_dayXYcoords.dat", function(y) -(y/500 * 180 - 90), function(x) (x/1000 * 360 - 180), TRUE)
}
# using splot to plot a 2D map from a data file in XY format, with the coordinates ready to use
Gprun("set pm3d map\nset size ratio -1;set lmargin 0;set rmargin 0;set tmargin 0;set bmargin 0;unset key;unset tics;unset border;unset colorbox\nset palette model RGB file \"" %s% testmap %s% 
    ".pal\" u ($1/255):($2/255):($3/255)\nset cbrange[0:" %s% NpaletteColors %s% "]\nset pm3d corners2color c1\nsplot \"earth_dayXYcoords.dat\" u 1:2:3 w pm3d notit", TRUE)

# plot the earth - 3d globe - matrix data with the indexed colors
Gprun("#set term pngcairo;set output \"globe6views.png\"\nset angles degrees\nunset colorbox\nunset border\nunset xtics\nunset ytics\nunset ztics\nset parametric\nset isosamples 13,13\nset mapping spherical\nset yzeroaxis linetype 0 linewidth 1.000 \nset urange [ -90.0000 : 90.0000 ] noreverse nowriteback\nset vrange [-180.0000 : 180.0000 ] noreverse nowriteback\nset hidden3d front offset 0 trianglepattern 3 undefined 1 altdiagonal bentover\nset view equal xyz\nset palette model RGB file \"" %s% 
    testmap %s% ".pal\" u ($1/255):($2/255):($3/255)\nset cbrange[0:" %s% NpaletteColors %s% "]\nset pm3d corners2color c1\nMwidth=" %s% Mwidth %s% "\nMheight=" %s% Mheight %s% "\nset multiplot layout 2,3\nset view 1,1,2,2 #North pole\nsplot \"" %s% 
    testmap %s% "matrixndx.dat\" matrix u (($1-Mheight)*180/Mheight):(-($2-Mheight/2)*180/Mheight):(1):3 w pm3d notit, cos(u)*cos(v),cos(u)*sin(v),sin(u) w l lc rgb \"grey\" notit, \"worldRmap.dat\" u 1:2:(1) w l lc rgb \"red\" notit\nset view 181,1,2,2 #South pole\nreplot\nset view 92,271,2,2 #Greenwich 180\nreplot\nset view 88,90,2,2 #Greenwich 0\nreplot\nset view 92,0,2,2 #Greenwich 270\nreplot\nset view 92,182,2,2 #Greenwich 90\nreplot\nunset multiplot", 
    TRUE)
GpImagePlot("globe6views.png")

# plot the earth - 3d globe - XY data with the indexed colors
Gprun("set angles degrees\nunset colorbox\nunset border\nunset xtics\nunset ytics\nunset ztics\nset parametric\nset mapping spherical\nset yzeroaxis linetype 0 linewidth 1.000 \nset vrange [ -90.0000 : 90.0000 ] noreverse nowriteback\nset urange [ -180.0000 : 180.0000 ] noreverse nowriteback\nset hidden3d front offset 0 trianglepattern 3 undefined 1 altdiagonal bentover\nset view equal xyz\nset isosamples 13,13\nset palette model RGB file \"" %s% 
    testmap %s% ".pal\" u ($1/255):($2/255):($3/255)\nset cbrange[0:" %s% NpaletteColors %s% "]\nset pm3d corners2color c1 \nMwidth=" %s% Mwidth %s% "\nMheight=" %s% Mheight %s% "-1\nset multiplot layout 2,3\nset view 358,2,2,2 #North pole\nsplot  \"" %s% 
    testmap %s% "XYndx.dat\" u (-(Mheight-$2)*180/Mheight):((Mheight/2-$1)*180/Mheight):(1):3 w pm3d notit, cos(u)*cos(v),cos(u)*sin(v),sin(u) w l lc rgb \"grey\" notit, \"worldRmap.dat\" u ($1):($2):(1) w l lc rgb \"red\" notit\nset view 178,178,2,2 #South pole\nreplot\nset view 88,188,2,2 #Greenwich 90\nreplot\nset view 92,90,2,2 #Greenwich 0\nreplot\nset view 88,358,2,2 #Greenwich 270\nreplot\nset view 88,268,2,2 #Greenwich 180\nreplot\nunset multiplot", 
    TRUE)


GpXYcoords2shpere(testmap %s% "XYndx.dat", testmap %s% "XYnewCOORDSndx.dat")
Gprun("set angles degrees\nunset colorbox\nunset border\nunset xtics\nunset ytics\nunset ztics\nset parametric\nset mapping spherical\nset yzeroaxis linetype 0 linewidth 1.000 \nset vrange [ -90.0000 : 90.0000 ] noreverse nowriteback\nset urange [ -180.0000 : 180.0000 ] noreverse nowriteback\nset hidden3d front offset 0 trianglepattern 3 undefined 1 altdiagonal bentover\nset view equal xyz\nset isosamples 13,13\nset palette model RGB file \"" %s% 
    testmap %s% ".pal\" u ($1/255):($2/255):($3/255)\nset cbrange[0:" %s% NpaletteColors %s% "]\nset pm3d corners2color c1\nset multiplot layout 2,3\nset view 2,2,2,2 #North pole\nsplot  \"" %s% 
    testmap %s% "XYnewCOORDSndx.dat\" u (-$2):($1):(1):3 w pm3d notit, cos(u)*cos(v),cos(u)*sin(v),sin(u) w l lc rgb \"grey\" notit, \"worldRmap.dat\" u 1:2:(1) w l lc rgb \"red\" notit\nset view 182,180,2,2 #South pole\nreplot\nset view 92,188,2,2 #Greenwich 90\nreplot\nset view 88,90,2,2 #Greenwich 0\nreplot\nset view 88,2,2,2 #Greenwich 270\nreplot\nset view 92,272,2,2 #Greenwich 180\nreplot\nunset multiplot", 
    TRUE)



 
