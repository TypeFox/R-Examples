## ----results = "hide"----------------------------------------------------
library(camtrapR)

## ------------------------------------------------------------------------

# load sample camera trap station table
data(camtraps)

# load sample record table
data(recordTableSample)

## ------------------------------------------------------------------------
Mapstest1 <- detectionMaps(CTtable     = camtraps,
                          recordTable  = recordTableSample,
                          Xcol         = "utm_x",
                          Ycol         = "utm_y",
                          stationCol   = "Station",
                          speciesCol   = "Species",
                          printLabels  = TRUE,
                          richnessPlot = TRUE,    # by setting this argument TRUE
                          speciesPlots = FALSE,
                          addLegend    = TRUE
)

## ------------------------------------------------------------------------
 # subset to 1 species
recordTableSample_PBE <- recordTableSample[recordTableSample$Species == "PBE",]
 
 Mapstest2 <- detectionMaps(CTtable      = camtraps,
                           recordTable   = recordTableSample_PBE,
                           Xcol          = "utm_x",
                           Ycol          = "utm_y",
                           stationCol    = "Station",
                           speciesCol    = "Species",
                           speciesToShow = "PBE",     # added
                           printLabels   = TRUE,
                           richnessPlot  = FALSE,     # changed
                           speciesPlots  = TRUE,      # changed
                           addLegend     = TRUE
)

## ----message=FALSE-------------------------------------------------------
# writing shapefiles requires packages rgdal and sp
 library(rgdal)
 library(sp)

## ------------------------------------------------------------------------
 # define shapefile name
 shapefileName       <- "recordShapefileTest"
 shapefileProjection <- "+proj=utm +zone=50 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
 
# run detectionMaps with shapefile creation
Mapstest3 <- detectionMaps(CTtable            = camtraps,
                          recordTable         = recordTableSample,
                          Xcol                = "utm_x",
                          Ycol                = "utm_y",
                          stationCol          = "Station",
                          speciesCol          = "Species",
                          richnessPlot        = FALSE,         # no richness plot
                          speciesPlots        = FALSE,         # no species plots
                          writeShapefile      = TRUE,          # but shaepfile creation
                          shapefileName       = shapefileName,
                          shapefileDirectory  = tempdir(),     # change this in your scripts!   
                          shapefileProjection = shapefileProjection
)
 
# check for the files that were created
list.files(tempdir(), pattern = shapefileName)

# load it as shapefile
shapefileTest <- readOGR(dsn = tempdir(), 
                         layer = shapefileName)

# we have a look at the attribute table
shapefileTest@data

# the output of detectionMaps is used as shapefile attribute table. Therefore, they are identical:
 all(shapefileTest@data == Mapstest3)
 

## ------------------------------------------------------------------------
# we first pick a species for our activity trials
species4activity <- "PBE"    # = Prionailurus bengalensis, Leopard Cat

## ------------------------------------------------------------------------

activityDensity(recordTable = recordTableSample,
                species     = species4activity)

## ------------------------------------------------------------------------
activityHistogram (recordTable = recordTableSample,
                   species     = species4activity)

## ------------------------------------------------------------------------
activityRadial(recordTable  = recordTableSample,
               species      = species4activity,
               lwd          = 3       # adjust line with of the plot
)

## ------------------------------------------------------------------------
activityRadial(recordTable       = recordTableSample,
               species           = species4activity,
               allSpecies        = FALSE,
               speciesCol        = "Species",
               recordDateTimeCol = "DateTimeOriginal",
               plotR             = TRUE,
               writePNG          = FALSE,
               lwd               = 3,
               rp.type           = "p"      # plot type = polygon
)

## ------------------------------------------------------------------------
# define species of interest
speciesA_for_activity <- "VTA"    # = Viverra tangalunga, Malay Civet
speciesB_for_activity <- "PBE"    # = Prionailurus bengalensis, Leopard Cat

# create activity overlap plot
activityOverlap (recordTable = recordTableSample,
                 speciesA    = speciesA_for_activity,
                 speciesB    = speciesB_for_activity,
                 writePNG    = FALSE,
                 plotR       = TRUE,
                 linecol     = c("red", "blue"),     # line colour of speciesA / speciesB
                 linewidth   = c(3,3),               # line width of speciesA / speciesB
                 add.rug     = TRUE
)


## ------------------------------------------------------------------------
reportTest <- surveyReport (recordTable          = recordTableSample,
                            CTtable              = camtraps,
                            speciesCol           = "Species",
                            stationCol           = "Station",
                            setupCol             = "Setup_date",
                            retrievalCol         = "Retrieval_date",
                            CTDateFormat         = "%d/%m/%Y", 
                            recordDateTimeCol    = "DateTimeOriginal",
                            recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                            CTHasProblems        = TRUE)

## ------------------------------------------------------------------------
str(reportTest)

## ------------------------------------------------------------------------
# here's the output of surveyReport

reportTest[[1]]    # camera trap operation times and image date ranges
reportTest[[2]]    # number of species by station
reportTest[[3]]    # number of events and number of stations by species
reportTest[[4]]    # number of species events by station
# reportTest[[5]] is identical to reportTest[[4]] except for the fact that it contains unobserved species with n_events = 0

