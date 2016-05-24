# 0. Load packages and specify working directory
library(RODBC); library(xlsx)
setwd("C:/aErer"); getwd(); dir()  # setting up global directory

# A. Read data in text format
daIns <- read.table(file = 'RawDataIns2.csv', header = TRUE, sep = ',')

# B. Read data in Excel format
# B1. Read one Excel sheet each time; output = a data frame
my.sheet <- c("dataImport", "dataExp", "nameImport", "nameExp", "source")
a1 <- read.xlsx(file = "RawDataAids.xlsx", sheetName = my.sheet[1])
a2 <- read.xlsx(file = "RawDataAids.xlsx", sheetName = my.sheet[2])

connectB <- odbcConnectExcel2007('RawDataAids.xlsx')  # open a connection
  sheet <- sqlTables(connectB); sheet$TABLE_NAME
  b1 <- sqlFetch(channel = connectB, sqtable = my.sheet[1])
  b2 <- sqlFetch(channel = connectB, sqtable = my.sheet[2])
odbcClose(connectB)  # close the connection

# B2. Read many sheets with a loop; output = a list
cc <- list()
for (k in my.sheet) {
  cc[[k]] <- read.xlsx(file = "RawDataAids.xlsx", sheetName = k)
}

dd <- list()
connectD <- odbcConnectExcel2007('RawDataAids.xlsx')
  sheet <- sqlTables(connectD); sheet$TABLE_NAME
  for (m in my.sheet) {
    dd[[m]] <- sqlFetch(channel = connectD, sqtable = m)
  }
odbcClose(connectD)

identical(a1, b1); identical(a2, b2); identical(cc[[3]], dd[[3]])
names(cc); cc[[1]][1:3, 1:7]

# C. Import graphics
library(png); library(jpeg)
nameA <- system.file("img", "Rlogo.png", package = "png" ); nameA
nameB <- system.file("img", "Rlogo.jpg", package = "jpeg"); nameB
# nameC <- "C:/bHome/EdPictures/MKSun2014Oct.jpg"  # your own picture

imageA <- readPNG(source  = nameA); dev.new(); grid.raster(imageA)
imageB <- readJPEG(source = nameB); dev.new(); grid.raster(imageB)
# imageC <- readJPEG(source = nameC); dev.new(); grid.raster(imageC)