### R code from vignette source 'rworldmapFAQ.Rnw'

###################################################
### code chunk number 1: rworldmapFAQ.Rnw:25-26
###################################################
options(width=70)


###################################################
### code chunk number 2: plotSetup
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")


###################################################
### code chunk number 3: rworldmapFAQ.Rnw:64-66
###################################################
require(rworldmap)
library(rworldmap)


###################################################
### code chunk number 4: showExampleCountryData
###################################################
data(countryExData)
countryExData[5:10,1:5]


###################################################
### code chunk number 5: joinCountryData2Map1
###################################################
data(countryExData)
sPDF <- joinCountryData2Map( countryExData,
                           , joinCode = "ISO3"
                           , nameJoinColumn = "ISO3V10" )


###################################################
### code chunk number 6: mapCountryData1
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
mapCountryData( sPDF, nameColumnToPlot="BIODIVERSITY" )


###################################################
### code chunk number 7: mapCountryData2
###################################################
mapParams <- mapCountryData( sPDF
                           , nameColumnToPlot="BIODIVERSITY"
                           , addLegend=FALSE )
do.call( addMapLegend, c(mapParams, legendWidth=0.5, legendMar = 2))


###################################################
### code chunk number 8: mapGriddedData1
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
data(gridExData)
mapGriddedData(gridExData)


###################################################
### code chunk number 9: mapHalfDegreeGridToCountries1
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
mapParams <- mapHalfDegreeGridToCountries(gridExData, addLegend=FALSE)
do.call( addMapLegend, c(mapParams, legendWidth=0.5, legendMar = 2))


###################################################
### code chunk number 10: country2Region1
###################################################
#Using country2Region to calculate mean Environmental Health index in Stern regions.
sternEnvHealth <- country2Region( inFile=countryExData
                                , nameDataColumn="ENVHEALTH"
                                , joinCode="ISO3"
                                , nameJoinColumn="ISO3V10"
                                , regionType="Stern"
                                , FUN="mean" )

print(sternEnvHealth)


###################################################
### code chunk number 11: mapByRegion1
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
mapByRegion( countryExData
           , nameDataColumn="CLIMATE"
           , joinCode="ISO3"
           , nameJoinColumn="ISO3V10"
           , regionType="Stern"
           , FUN="mean" )


###################################################
### code chunk number 12: colourPalette
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")

#joining the data to a map
sPDF <- joinCountryData2Map( countryExData
                           , joinCode = "ISO3"
                           , nameJoinColumn = "ISO3V10"
                           )
                                       
#creating a user defined colour palette
op <- palette(c('green','yellow','orange','red')) 
#find quartile breaks
cutVector <- quantile(sPDF@data[["BIODIVERSITY"]],na.rm=TRUE)
#classify the data to a factor
sPDF@data[["BIOcategories"]] <- cut( sPDF@data[["BIODIVERSITY"]]
                                   , cutVector
                                   , include.lowest=TRUE )
#rename the categories
levels(sPDF@data[["BIOcategories"]]) <- c('low', 'med', 'high', 'vhigh')
#mapping
mapCountryData( sPDF
              , nameColumnToPlot='BIOcategories'
              , catMethod='categorical'
              , mapTitle='Biodiversity categories'
              , colourPalette='palette'
              , oceanCol='lightblue'
              , missingCountryCol='white' )


###################################################
### code chunk number 13: finalFigure2
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
         
mapCountryData( sPDF
              , nameColumnToPlot='BIOcategories'
              , catMethod='categorical'
              , mapTitle='Biodiversity categories'
              , colourPalette='palette'
              , oceanCol='lightblue'
              , missingCountryCol='white' 
              , mapRegion='Eurasia'
              , borderCol='black' )
  



###################################################
### code chunk number 14: selectedCountries
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
sPDF <- getMap()  
#select countries from the map       
sPDF <-sPDF[which(sPDF$LDC=='LDC'),] 
mapCountryData( sPDF
              , nameColumnToPlot='continent'
              , colourPalette='rainbow'
              , mapTitle='Least Developed Countries' )




###################################################
### code chunk number 15: bubblePlot
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")

mapBubbles( dF=getMap()
          , nameZSize="POP_EST"
          , nameZColour="continent"
          , colourPalette='rainbow'
          , oceanCol='lightblue'
          , landCol='wheat' ) 


###################################################
### code chunk number 16: mtext
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")

sPDF <- getMap()  
#select countries from the map       
sPDF <-sPDF[which(sPDF$continent=='Africa'),] 
mapBubbles( dF=getMap()
          , nameZSize="POP_EST"
          , nameZColour="continent"
					, mapTitle='Population'
					, addColourLegend = FALSE) 
          
mtext("Source: Andy South, The R Journal Vol. 3/1, June 2011",side=1,line=-1)          
          


###################################################
### code chunk number 17: projection
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")

library(rgdal)
#first get countries excluding Antarctica which crashes spTransform
sPDF <- getMap()[-which(getMap()$ADMIN=='Antarctica'),]
#transform to robin for the Robinson projection
sPDF <- spTransform(sPDF, CRS=CRS("+proj=robin +ellps=WGS84")) 
       
mapCountryData( sPDF
          , nameColumnToPlot="REGION"
					, mapTitle='Robinson Projection'
					, colourPalette='topo'
					, addLegend = FALSE)                    
          


###################################################
### code chunk number 18: classInt_RColorBrewer
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")

library(classInt)
library(RColorBrewer)

#getting example data and joining to a map
data("countryExData",envir=environment(),package="rworldmap")
sPDF <- joinCountryData2Map( countryExData
                           , joinCode = "ISO3"
                           , nameJoinColumn = "ISO3V10"
                           , mapResolution='coarse' )

#getting class intervals using a 'jenks' classification in classInt package
classInt <- classIntervals( sPDF[["EPI"]], n=5, style="jenks")
catMethod = classInt[["brks"]]

#getting a colour scheme from the RColorBrewer package
colourPalette <- brewer.pal(5,'RdPu')

#calling mapCountryData with the parameters from classInt and RColorBrewer
mapParams <- mapCountryData( sPDF
                           , nameColumnToPlot="EPI"
                           , addLegend=FALSE
                           , catMethod = catMethod
                           , colourPalette=colourPalette )
do.call( addMapLegend
       , c( mapParams
          , legendLabels="all"
          , legendWidth=0.5
          , legendIntervals="data"
          , legendMar = 2 ))


###################################################
### code chunk number 19: margins
###################################################
oldPar <- par(mar=c(0, 0, 0, 0)) 
par(oldPar)


###################################################
### code chunk number 20: layout1
###################################################

#set margins to zero for the subplots
oldPar <- par(mar=c(0, 0, 0, 0)) 
nPanels <- layout( cbind(c(0,1,2,3,4,5),c(0,6,7,8,9,10))
                 , heights=c(lcm(0.5),1,1,1,1,1)
                 , respect=F)      
       
layout.show(nPanels)
par(oldPar)


###################################################
### code chunk number 21: layoutMonthly
###################################################

#set margins to zero for the subplots
oldPar <- par(mar=c(0, 0, 0, 0)) 
nPanels <- layout( rbind(c(0,0,0),c(1,2,3),c(4,5,6),c(7,8,9),c(10,11,12))
                 , heights=c(lcm(0.5),1,1,1,1)
                 , respect=F )    
       
layout.show(nPanels)
par(oldPar)


###################################################
### code chunk number 22: linesLatLon
###################################################
abline(h=0) 
abline(v=0)   
abline(h=c(-20,20),lty=2,col='grey')  


