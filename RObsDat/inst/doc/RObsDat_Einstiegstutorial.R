### R code from vignette source 'RObsDat_Einstiegstutorial.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: default (eval = FALSE)
###################################################
## getDefaultDB()


###################################################
### code chunk number 2: sqlite (eval = FALSE)
###################################################
## require("RObsDat")
## require("RSQLite")
## m <- dbDriver("SQLite")
## dbname = "database.db"
## con <- dbConnect(m, dbname = dbname)
## sqhandler <-  new("odm1_1Ver", con=con)
## options(odm.handler=sqhandler)


###################################################
### code chunk number 3: engines (eval = FALSE)
###################################################
## #connect to postgreSQL database
## require("RObsDat")
## require("RPostgreSQL")
## m <- dbDriver("PostgreSQL")
## con <- dbConnect(m, user="a_user", password="secret", dbname="obsdat")
## sqhandler <-  new("odm1_1Ver", con=con)
## options(odm.handler=sqhandler)
## 
## #connect to MySQL database
## require("RObsDat")
## require("RMySQL")
## m <- dbDriver("MySQL")
## con <- dbConnect(m, user="a_user", password="secret", dbname="obsdat")
## sqhandler <-  new("odm1_1Ver", con=con)
## options(odm.handler=sqhandler)


###################################################
### code chunk number 4: part1
###################################################
require("RObsDat")
getDefaultDB()


###################################################
### code chunk number 5: part2
###################################################
#Store metadata in database
getMetadata("SpatialReference", SRSName="WGS84", exact=TRUE)

addSite(Code="testSpatialPoints", Name="Virtual test site", x=25, y=56,
	LatLongDatum="WGS84", Elevation=350, State="Germany")

getMetadata("Units", Name="degree celsius")
getMetadata("VariableName", Term="Temperature")

addVariable(Name="Temperature, transducer signal", Unit="degree celsius", ValueType="Field Observation",
	GeneralCategory="Hydrology", Code="test_temp")

addQualityControlLevel(ID=2,Code="test_ok", Definition="The default values")

addISOMetadata(TopicCategory="Unknown", Title="Testdata",
	Abstract="This data is created to test the functions of RObsDat")
addSource(Organization="Your Org", SourceDescription="Madeup data", 
	SourceLink="RObsDat Documentation", ContactName="Yourself",
	Metadata="Testdata")


###################################################
### code chunk number 6: part2
###################################################
require(xts)

example.data <- xts(1:40, seq(as.POSIXct("2014-01-01", tz="UTC"), 
		as.POSIXct("2014-02-09", tz="UTC"), length.out=40))
example.data[40] <- 30
example.data[35] <- 22


###################################################
### code chunk number 7: part3
###################################################
addDataValues(example.data[1:20], Site="Virtual test site", Variable="test_temp",  
	Source="Madeup", QualityControlLevel="test_ok")



###################################################
### code chunk number 8: part3
###################################################
#Avoid duplicates autmatically
example.data[15] <- 30
addDataValues(example.data, Site="Virtual test site", Variable="test_temp",  
	Source="Madeup", QualityControlLevel="test_ok")


###################################################
### code chunk number 9: part6
###################################################
allData <- getDataValues()
testSiteData <- getDataValues(Site="test")


###################################################
### code chunk number 10: part5
###################################################
#Version management

testSiteData <- getDataValues(Site="test")

to.correct <- which(testSiteData@data > 30)
testSiteData@data[to.correct,] <- 20
testSiteData@data[39,] <- 32

#ToDo - doesn't pass test
#updateDataValues(testSiteData, "Correction of wrong value")

ver2 <- testSiteData
ver2@data[10:13,] <- 60
#updateDataValues(ver2, "Changing more data")

testSiteData <- getDataValues(Site="test")
ver3 <- testSiteData
ver3@data[30:32,] <- 33
#updateDataValues(ver3, "Ups, I used 60 instead of 33 by mistake")



###################################################
### code chunk number 11: partxyz
###################################################
getDataVersions()

versionQuery1 <- getDataValues(Site=1, VersionID=1)
#stplot(versionQuery1, mode="ts")

versionQuery2 <- getDataValues(Site=1, VersionID=2)
#stplot(versionQuery2, mode="ts")


