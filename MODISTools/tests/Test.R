# Testing for internet connectivity, the connection to the MODIS SOAP WSDL Server and it's Web Service
# Description Language, for the XML response from the Web Service method, and for the functions of
# MODISTools.

# Load data to be used for testing.
rm(list = ls())
library(MODISTools)
data(SubsetExample, FindIDExample, QualityCheckExample, TransectExample, EndCoordinatesExample, ConvertExample)
library(RCurl)  # Will use some RCurl and XML functions explicitly in testing.
library(XML)

options(warn = 2)

## Following lines of code testing for internet connectivity and server access, are from
## R testing: .../tests/internet.R
# Check for internet capability.
if(!capabilities("http/ftp")) q()

# Check for internet connectivity.
if(.Platform$OS.type == "unix" && is.null(nsl("cran.r-project.org"))) q()

# Check we can reach the server for lpdaac modis web service.
if(.Platform$OS.type == "unix" && is.null(nsl("daac.ornl.gov"))) q()

# Check the web service is currently responsive.
if(class(try(GetProducts(), silent = TRUE)) == "try-error") q()
##

# Check the XML response is as expected.
getsubset.xml <- paste('
<soapenv:Envelope xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xsd="http://www.w3.org/2001/XMLSchema"
              xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/" xmlns:mod="http://daac.ornl.gov/MODIS_webservice">
                        <soapenv:Header/>
                        <soapenv:Body>
                        <mod:getsubset soapenv:encodingStyle="http://schemas.xmlsoap.org/soap/encoding/">
                        <Latitude xsi:type="xsd:float">', 51.41363, '</Latitude>
                        <Longitude xsi:type="xsd:float">', -0.64875, '</Longitude>
                        <Product xsi:type="xsd:string">', "MOD13Q1", '</Product>
                        <Band xsi:type="xsd:string">', "250m_16_days_EVI", '</Band>
                        <MODIS_Subset_Start_Date xsi:type="xsd:string">', "A2001001", '</MODIS_Subset_Start_Date>
                        <MODIS_Subset_End_Date xsi:type="xsd:string">', "A2001025", '</MODIS_Subset_End_Date>
                        <Km_Above_Below xsi:type="xsd:string">', 0, '</Km_Above_Below>
                        <Km_Left_Right xsi:type="xsd:string">', 0, '</Km_Left_Right>
                        </mod:getsubset>
                        </soapenv:Body>
                        </soapenv:Envelope>',
                        sep = "")

header.fields <- c(Accept = "text/xml",
                   Accept = "multipart/*",
                   'Content-Type' = "text/xml; charset=utf-8",
                   SOAPAction = "")

reader <- basicTextGatherer()
header <- basicTextGatherer()

curlPerform(url = "http://daac.ornl.gov/cgi-bin/MODIS/GLBVIZ_1_Glb_subset/MODIS_webservice.pl",
            httpheader = header.fields,
            postfields = getsubset.xml,
            writefunction = reader$update,
            verbose = FALSE)

# Check the server is not down by insepcting the XML response for internal server error message.
if(grepl("Internal Server Error", reader$value())) q()

xmlRoot(xmlTreeParse(reader$value()))
###

# Check FindID example
FindID(ID = SubsetExample, Data = FindIDExample)

# Check QualityCheck example
EVIdata <- QualityCheckExample[1:5, ]
QAdata <- QualityCheckExample[6:10, ]

QualityCheck(Data = EVIdata, Product = "MOD13Q1", Band = "250m_16_days_EVI", NoDataFill = -3000,
             QualityBand = "250m_16_days_pixel_reliability", QualityScores = QAdata, QualityThreshold = 0)
###

# Check MODIS subset uses this output to produce correctly downloaded files.
request <- GetSubset(Lat = SubsetExample$lat, Long = SubsetExample$long, Product = "MCD12Q1", Band = "Land_Cover_Type_1",
                     StartDate = "A2005001", EndDate = "A2005001", KmAboveBelow = 0, KmLeftRight = 0)$subset[1]
if(grepl("Server is busy handling other requests", request) | grepl("System overloaded", request) |
     grepl("Downloading from the web service is currently not working", request)){
  q()
} else {
  # Check GetSubset is producing the correct output.
  # Use GetProducts, GetBands, and GetDates, to specify the GetSubset request.
  Product <- GetProducts()[1]
  Band <- GetBands(Product)[1]
  Dates <- GetDates(SubsetExample$lat, SubsetExample$long, Product)[1:2]

  GetSubset(Lat = SubsetExample$lat, Long = SubsetExample$long, Product = Product, Band = Band,
            StartDate = Dates[1], EndDate = Dates[1], KmAboveBelow = 0, KmLeftRight = 0)

  MODISSubsets(LoadDat = SubsetExample, Product = "MCD12Q1", Bands = c("Land_Cover_Type_1"),
               Size = c(1,1), StartDate = TRUE)

  MODISSummaries(LoadDat = SubsetExample, Product = "MCD12Q1", Band = "Land_Cover_Type_1",
                 ValidRange = c(0,254), NoDataFill = 255, ScaleFactor = 1, StartDate = TRUE)
}

# Check the MODISSummaries file outputs are consistent.
SummaryFile <- read.csv(list.files(pattern = "MODIS_Summary"))
DataFile <- read.csv(list.files(pattern = "MODIS_Data"))
file.check <- all(SummaryFile$mean.band == DataFile[1,which(grepl("pixel", names(DataFile)))])
if(is.na(file.check)){
  warning("The two output files from MODISSummaries are not consistent.")
}
if(!file.check){
  warning("The two output files from MODISSummaries are not consistent.")
}

# Check example of MODISTransects
request <- GetSubset(Lat = SubsetExample$lat, Long = SubsetExample$long, Product = "MOD13Q1", Band = "250m_16_days_EVI",
                     StartDate = "A2000049", EndDate = "A2000049", KmAboveBelow = 0, KmLeftRight = 0)$subset[1]
if(grepl("Server is busy handling other requests", request) | grepl("System overloaded", request) |
     grepl("Downloading from the web service is currently not working", request)){
  q()
} else {
  MODISTransects(LoadData = TransectExample, Product = "MCD12Q1", Bands = c("Land_Cover_Type_1"),
                 Size = c(0,0), StartDate = TRUE)
}

# Check EndCoordinates example
EndCoordinates(LoadDat = EndCoordinatesExample, Distance = 2000, Angle = 90, AngleUnits = "degrees")

# Check ConvertToDD example
ConvertToDD(XY = ConvertExample, LatColName = "lat", LongColName = "long")

# Check ExtractTile example
TileExample <- read.csv(list.files(pattern = "MODIS_Data"))
TileExample <- TileExample[ ,which(grepl("pixel", names(TileExample)))]

dim(TileExample)
dim(ExtractTile(Data = TileExample, Rows = c(5,1), Cols = c(5,1), Grid = FALSE))
ExtractTile(Data = TileExample, Rows = c(5,1), Cols = c(5,1), Grid = FALSE)

matrix(TileExample, nrow = 5, ncol = 5, byrow = TRUE)
ExtractTile(Data = TileExample, Rows = c(5,1), Cols = c(5,1), Grid = TRUE)

# Check LandCover on previously downloaded data from MODISSubsets
LandCover(Band = "Land_Cover_Type_1")

rm(list = ls())
options(warn = 0)