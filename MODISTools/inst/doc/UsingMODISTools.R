### R code from vignette source 'UsingMODISTools.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: UsingMODISTools.Rnw:25-36
###################################################
library(MODISTools)

# Makes copy-paste much less painful
options(continue = ' ')
options(width = 90)
options(prompt = '> ')

options(SweaveHooks = list(fig=function() par(mgp=c(2.5,1,0),
                                              mar=c(4,4,2,1),
                                              oma=c(0,0,1,0),
                                              cex.main=0.8)))


###################################################
### code chunk number 2: UsingMODISTools.Rnw:46-48
###################################################
data(ConvertExample)
ConvertExample


###################################################
### code chunk number 3: UsingMODISTools.Rnw:51-55
###################################################
modis.subset <-
  ConvertToDD(XY = ConvertExample, LatColName = "lat", LongColName = "long")
modis.subset <- data.frame(lat = modis.subset[ ,1], long = modis.subset[ ,2])
modis.subset


###################################################
### code chunk number 4: UsingMODISTools.Rnw:58-60
###################################################
modis.subset$start.date <- rep(2003, nrow(modis.subset))
modis.subset$end.date <- rep(2006, nrow(modis.subset))


###################################################
### code chunk number 5: UsingMODISTools.Rnw:67-68 (eval = FALSE)
###################################################
## GetProducts()


###################################################
### code chunk number 6: UsingMODISTools.Rnw:70-73
###################################################
c("MCD12Q1", "MCD12Q2", "MCD43A1", "MCD43A2", "MCD43A4", "MOD09A1",
"MOD11A2", "MOD13Q1", "MOD15A2", "MOD15A2GFS", "MOD16A2", "MOD17A2_51",
"MOD17A3", "MYD09A1", "MYD11A2", "MYD13Q1", "MYD15A2")


###################################################
### code chunk number 7: UsingMODISTools.Rnw:75-76 (eval = FALSE)
###################################################
## GetBands(Product = "MOD13Q1")


###################################################
### code chunk number 8: UsingMODISTools.Rnw:78-84
###################################################
c("250m_16_days_blue_reflectance", "250m_16_days_MIR_reflectance",
"250m_16_days_NIR_reflectance", "250m_16_days_pixel_reliability",
"250m_16_days_red_reflectance", "250m_16_days_relative_azimuth_angle",
"250m_16_days_sun_zenith_angle", "250m_16_days_view_zenith_angle",
"250m_16_days_VI_Quality", "250m_16_days_NDVI",
"250m_16_days_EVI", "250m_16_days_composite_day_of_the_year")


###################################################
### code chunk number 9: UsingMODISTools.Rnw:89-90 (eval = FALSE)
###################################################
## GetDates(Product = "MOD13Q1", Lat = modis.subset$lat[1], Long = modis.subset$long[1])


###################################################
### code chunk number 10: UsingMODISTools.Rnw:96-99 (eval = FALSE)
###################################################
## MODISSubsets(LoadDat = modis.subset, Products = "MOD13Q1",
##              Bands = c("250m_16_days_EVI", "250m_16_days_pixel_reliability"),
##              Size = c(1,1))


###################################################
### code chunk number 11: UsingMODISTools.Rnw:104-107 (eval = FALSE)
###################################################
## subset.string <- read.csv(list.files(pattern = ".asc")[1],
##                                 header = FALSE, as.is = TRUE)
## subset.string[1, ]


###################################################
### code chunk number 12: UsingMODISTools.Rnw:109-113
###################################################
subset.string <- read.csv(paste("./MODISSubsetsMOD13Q1/",
                                list.files(path = "./MODISSubsetsMOD13Q1", pattern = ".asc")[1]
                                , sep = ""), header = FALSE, as.is = TRUE)
subset.string[1, ]


###################################################
### code chunk number 13: UsingMODISTools.Rnw:123-127 (eval = FALSE)
###################################################
## MODISSummaries(LoadDat = modis.subset, Product = "MOD13Q1", Bands = "250m_16_days_EVI",
##                ValidRange = c(-2000,10000), NoDataFill = -3000, ScaleFactor = 0.0001,
##                QualityScreen = TRUE, QualityBand = "250m_16_days_pixel_reliability",
##                QualityThreshold = 0)


###################################################
### code chunk number 14: UsingMODISTools.Rnw:133-135 (eval = FALSE)
###################################################
## TileExample <- read.csv(list.files(pattern = "MODIS_Data"))
## TileExample <- TileExample[ ,which(grepl("250m_16_days_EVI", names(TileExample)))]


###################################################
### code chunk number 15: UsingMODISTools.Rnw:137-141
###################################################
TileExample <- read.csv(paste("./MODISSummaries/",
                              list.files(path = "./MODISSummaries/",
                                         pattern = "Data"), sep = ""))
TileExample <- TileExample[ ,which(grepl("250m_16_days_EVI", names(TileExample)))]


###################################################
### code chunk number 16: UsingMODISTools.Rnw:144-148
###################################################
dim(TileExample)
dim(ExtractTile(Data = TileExample, Rows = c(9,2), Cols = c(9,2), Grid = FALSE))
head(ExtractTile(Data = TileExample, Rows = c(9,2), Cols = c(9,2), Grid = FALSE),
     n = 2)


###################################################
### code chunk number 17: UsingMODISTools.Rnw:151-153
###################################################
matrix(TileExample[1, ], nrow = 9, ncol = 9, byrow = TRUE)
ExtractTile(Data = TileExample, Rows = c(9,2), Cols = c(9,2), Grid = TRUE)[ , ,1]


###################################################
### code chunk number 18: UsingMODISTools.Rnw:159-163 (eval = FALSE)
###################################################
## dir.create('./LandCover')
## setwd('./LandCover')
## MODISSubsets(LoadDat = modis.subset, Product = "MCD12Q1", Bands = "Land_Cover_Type_1",
##              Size = c(1,1))


###################################################
### code chunk number 19: UsingMODISTools.Rnw:166-170 (eval = FALSE)
###################################################
## LandCover(Band = "Land_Cover_Type_1")
## 
## land.summary <- read.csv(list.files(pattern = "MODIS_Land_Cover_Summary"))
## head(land.summary)


###################################################
### code chunk number 20: UsingMODISTools.Rnw:172-177
###################################################
land.summary <- read.csv(paste("./LandCover/",
                               list.files(path = "./LandCover/",
                                          pattern = "LandCoverSummary"),
                               sep = ""))
head(land.summary)


