library(maps)
library(mapdata)
# first make sure that all the data files exist
if ((!file.exists("NOAACoastline.dat")) | (!file.exists("worldmer.dat")) | (!file.exists("worldpar.dat")) | (!file.exists("tissot.dat")) | (!file.exists("earth_dayXYcoords.dat")) | 
    (!file.exists("earth_day.pal"))) {
    stop("Please run demo(maps2) before this demo, so that all the data files are created")
}

# now show the examples
vmaps <- c("NOAACoastline.dat", "worldmer.dat", "worldpar.dat")
GpPlotMap(vmaps, linestyle = c(1, 2, 3))
GpPlotMap(vmaps, linestyle = c(1, 2, 3), projection = "Mercator")

GpPlotMap(projection = "Mercator", maprastfile = "earth_dayXYcoords.dat", maprastpalette = "earth_day.pal")
GpPlotMap(c("worldmer.dat", "worldpar.dat"), projection = "Mercator", linestyle = c(1, 2, 3), maprastfile = "earth_dayXYcoords.dat", maprastpalette = "earth_day.pal")

# create meridian and parallel lines for plot
GpMapMerpar("worldpar15.dat", "worldmer15.dat", 15, 15)
GpMapMerpar("worldpar20.dat", "worldmer20.dat", 20, 20)

vmaps2 <- c("NOAACoastline.dat", "worldmer15.dat", "worldpar15.dat", "tissot.dat")
vstyle2 <- c(1, 2, 3, 5)
GpPlotMap(vmaps2, linestyle = vstyle2, projection = "Mercator")  # ,AdditionalCode='set term png;set output 'Mercator1.png'\n'
GpPlotMap(vmaps2, linestyle = vstyle2, projection = "Ortographic")  # ,AdditionalCode='set term png;set output 'Ortographic1.png'\n'
GpPlotMap(vmaps2, linestyle = vstyle2, projection = "EstereoAzimutal")  # ,AdditionalCode='set term png;set output 'EstereoAzimutal1.png'\n'
GpPlotMap(vmaps2, linestyle = vstyle2, projection = "PlateCarree")  # ,AdditionalCode='set term png;set output 'PlateCarree1.png'\n'
GpPlotMap(vmaps2, linestyle = vstyle2, projection = "Lambert")  # ,AdditionalCode='set term png;set output 'Lambert1.png'\n'
GpPlotMap(vmaps2, linestyle = vstyle2, projection = "SansonFlamsteed")  # ,AdditionalCode='set term png;set output 'SansonFlamsteed1.png'\n'
GpPlotMap(vmaps2, linestyle = vstyle2, projection = "AlbersConical")  # ,AdditionalCode='set term png;set output 'AlbersConical1.png'\n'
GpPlotMap(vmaps2, linestyle = vstyle2, projection = "EckertI")  # ,AdditionalCode='set term png;set output 'EckertI1.png'\n'
GpPlotMap(vmaps2, linestyle = vstyle2, projection = "HammerWagner")  # ,AdditionalCode='set term png;set output 'HammerWagner1.png'\n'
GpPlotMap(vmaps2, linestyle = vstyle2, projection = "WernersEquivalent")  # ,AdditionalCode='set term png;set output 'WernersEquivalent1.png'\n'
GpPlotMap(vmaps2, linestyle = vstyle2, projection = "NaturalEarth")  # ,AdditionalCode='set term png;set output 'NaturalEarth1.png'\n'
GpPlotMap(vmaps2, linestyle = vstyle2, projection = "Aitoff")  # ,AdditionalCode='set term png;set output 'Aitoff1.png'\n'
GpPlotMap(vmaps2, linestyle = vstyle2, projection = "Winkeltripel")  # ,AdditionalCode='set term png;set output 'Winkeltripel1.png'\n'
GpPlotMap(vmaps2, linestyle = vstyle2, projection = "Robinson")  # ,AdditionalCode='set term png;set output 'Robinson1.png'\n'

vmaps3 <- c("worldmer15.dat", "worldpar15.dat")
vstyle3 <- c(2, 3)
GpPlotMap(vmaps3, projection = "Mercator", linestyle = vstyle3, maprastfile = "earth_dayXYcoords.dat", maprastpalette = "earth_day.pal")  # ,AdditionalCode='set term pngcairo;set output 'Mercator2.png'\n'
GpPlotMap(vmaps3, projection = "Ortographic", linestyle = vstyle3, maprastfile = "earth_dayXYcoords.dat", maprastpalette = "earth_day.pal")  # ,AdditionalCode='set term pngcairo;set output 'Ortographic2.png'\n'
GpPlotMap(vmaps3, projection = "EstereoAzimutal", linestyle = vstyle3, maprastfile = "earth_dayXYcoords.dat", maprastpalette = "earth_day.pal")  # ,AdditionalCode='set term pngcairo;set output 'EstereoAzimutal2.png'\n'
GpPlotMap(vmaps3, projection = "PlateCarree", linestyle = vstyle3, maprastfile = "earth_dayXYcoords.dat", maprastpalette = "earth_day.pal")  # ,AdditionalCode='set term pngcairo;set output 'PlateCarree2.png'\n'
GpPlotMap(vmaps3, projection = "Lambert", linestyle = vstyle3, maprastfile = "earth_dayXYcoords.dat", maprastpalette = "earth_day.pal")  # ,AdditionalCode='set term pngcairo;set output 'Lambert2.png'\n'
GpPlotMap(vmaps3, projection = "SansonFlamsteed", linestyle = vstyle3, maprastfile = "earth_dayXYcoords.dat", maprastpalette = "earth_day.pal")  # ,AdditionalCode='set term pngcairo;set output 'SansonFlamsteed2.png'\n'
GpPlotMap(vmaps3, projection = "AlbersConical", linestyle = vstyle3, maprastfile = "earth_dayXYcoords.dat", maprastpalette = "earth_day.pal")  # ,AdditionalCode='set term pngcairo;set output 'AlbersConical2.png'\n'
GpPlotMap(vmaps3, projection = "EckertI", linestyle = vstyle3, maprastfile = "earth_dayXYcoords.dat", maprastpalette = "earth_day.pal")  # ,AdditionalCode='set term pngcairo;set output 'EckertI2.png'\n'
GpPlotMap(vmaps3, projection = "HammerWagner", linestyle = vstyle3, maprastfile = "earth_dayXYcoords.dat", maprastpalette = "earth_day.pal")  # ,AdditionalCode='set term pngcairo;set output 'HammerWagner2.png'\n'
GpPlotMap(vmaps3, projection = "WernersEquivalent", linestyle = vstyle3, maprastfile = "earth_dayXYcoords.dat", maprastpalette = "earth_day.pal")  #,AdditionalCode='set term pngcairo;set output 'WernersEquivalent2.png'\n' 
GpPlotMap(vmaps3, projection = "NaturalEarth", linestyle = vstyle3, maprastfile = "earth_dayXYcoords.dat", maprastpalette = "earth_day.pal")  # ,AdditionalCode='set term pngcairo;set output 'NaturalEarth2.png'\n'
GpPlotMap(vmaps3, projection = "Aitoff", linestyle = vstyle3, maprastfile = "earth_dayXYcoords.dat", maprastpalette = "earth_day.pal")  # ,AdditionalCode='set term pngcairo;set output 'Aitoff2.png'\n'
GpPlotMap(vmaps3, projection = "Robinson", linestyle = vstyle3, maprastfile = "earth_dayXYcoords.dat", maprastpalette = "earth_day.pal")  # ,AdditionalCode='set term pngcairo;set output 'Robinson2.png'\n'
GpPlotMap(vmaps3, projection = "Winkeltripel", linestyle = vstyle3, maprastfile = "earth_dayXYcoords.dat", maprastpalette = "earth_day.pal")  # ,AdditionalCode='set term pngcairo;set output 'Winkeltripel2.png'\n'


# Ortographic Projection tangent plane at latitude 90 N and at latitude 90 S
vmaps2 <- c("worldmer15.dat", "worldpar15.dat", "tissot.dat", "NOAACoastline.dat")
vstyle2 <- c(2, 3, 5, 1)
Gprun("#set term png;set output \"NOAACoastline_multi1.png\"\nload \"projections.gnu\"\nset multiplot layout 1,2 scale 1.5 title \"Ortographic Projection\"\n" %s% GpPlotMap(vmaps2, 
    linestyle = vstyle2, projection = "Ortographic", plotTitle = "", AdditionalCode = "set title \"Northern Hemisphere\" offset 0,-2", projectionInit = "90,0,90,180", returnCode = TRUE) %s% 
    GpPlotMap(vmaps2, linestyle = vstyle2, projection = "Ortographic", plotTitle = "Southern Hemisphere", projectionInit = "-90,0,90,180", returnCode = TRUE) %s% "\nunset multiplot", 
    1)

# Ortographic Projection tangent plane at the Equator and Greenwich, plus a scaled view of Northern Europe
Gprun("#set term png;set output \"NOAACoastline_multi2.png\"\nload \"projections.gnu\"\nset multiplot layout 1,2 scale 1.5 title \"Ortographic Projection\"\n" %s% GpPlotMap(vmaps2, 
    linestyle = vstyle2, projection = "Ortographic", plotTitle = "", AdditionalCode = "set title \"Equator and Greenwich\" offset 0,-2", projectionInit = "0,0,90,90", returnCode = TRUE) %s% 
    GpPlotMap(c("worldmer.dat", "worldpar.dat", "NOAACoastline.dat"), linestyle = c(2, 3, 1), projection = "Ortographic", plotTitle = "Northern Europe", AdditionalCode = "ap=100", projectionInit = "60,20,10,20", 
        returnCode = TRUE) %s% "\nunset multiplot", 1)

# Albers Conical Equivalent Projection showing the difference between centering at northern versus southern latitudes.
Gprun("#set term png;set output \"NOAACoastline_multi3.png\"\nload \"projections.gnu\"\nset multiplot layout 1,2 title \"Albers Conical Equivalent Projection\"\n" %s% GpPlotMap(vmaps2, 
    linestyle = vstyle2, projection = "AlbersConical", plotTitle = "", AdditionalCode = "set title \"latitude 20 S\" offset 0,-2", projectionInit = "-20,0", returnCode = TRUE) %s% GpPlotMap(vmaps2, 
    linestyle = vstyle2, projection = "AlbersConical", plotTitle = "latitude 20 N", AdditionalCode = "lat0=20", projectionInit = "20,0", returnCode = TRUE) %s% "\nunset multiplot", 
    1)

 
