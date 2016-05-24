### R code from vignette source 'rworldmap.rnw'

###################################################
### code chunk number 1: rworldmap.rnw:70-71
###################################################
library(rworldmap)


###################################################
### code chunk number 2: rworldmap.rnw:82-83
###################################################
options(width=70)


###################################################
### code chunk number 3: plotSetup
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")


###################################################
### code chunk number 4: eg1
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")

data("countryExData",envir=environment())
sPDF <- joinCountryData2Map(countryExData,
              joinCode = "ISO3",
              nameJoinColumn = "ISO3V10",
              mapResolution = "coarse"
              )
#map the data with no legend
mapParams <- mapCountryData( sPDF,
              nameColumnToPlot="BIODIVERSITY",
              addLegend=FALSE,
              oceanCol="lightblue",
              missingCountryCol="grey"
              )
#add a modified legend
#do.call( addMapLegend, c(mapParams,
#            legendLabels="all",
#            legendWidth=0.5,
#            legendMar = 2
#            ))

#colourVector = mapParams[["colourVector"]]
#cutVector = mapParams[["cutVector"]]
#warning(paste("length colourVector=",length(mapParams$colourVector),"\n"
#             ,"length cutVector=",length(mapParams$cutVector),"\n"
#						 ,"length plottedData=",length(mapParams$plottedData),"\n"
#						 ,"catMethod=",mapParams$catMethod,"\n"
#						 ,"colourPalette=",mapParams$colourPalette,"\n"
#						 ))						 
#warning(paste("length colourVector=",length(colourVector),"length cutVector=",length(cutVector),"\n"))

#6/11/12 problem here
addMapLegend(colourVector = mapParams[["colourVector"]],
             cutVector = mapParams[["cutVector"]], 
             legendLabels="all",
             legendWidth=0.5,
             legendMar = 2
             )



###################################################
### code chunk number 5: eg2
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
data("gridExData",envir=environment())
mapParams <- mapGriddedData(oceanCol="lightblue",borderCol=NA,landCol="white",addLegend=FALSE, colourPalette="topo")
#do.call( addMapLegend, c(mapParams,legendLabels="all",legendWidth=0.5, legendMar = 2))
addMapLegend(colourVector = mapParams[["colourVector"]],
             cutVector = mapParams[["cutVector"]],   
             legendLabels="all",
             legendWidth=0.5,
             legendMar = 2
             )

title("Population per half degree grid cell")


###################################################
### code chunk number 6: showExampleCountryData
###################################################
data(countryExData)
countryExData[5:10,1:5]


###################################################
### code chunk number 7: joinCountryData2Map1
###################################################
data(countryExData)
sPDF <- joinCountryData2Map( countryExData
                           , joinCode = "ISO3"
                           , nameJoinColumn = "ISO3V10")


###################################################
### code chunk number 8: mapCountryData1
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
mapCountryData( sPDF, nameColumnToPlot="BIODIVERSITY" )


###################################################
### code chunk number 9: mapCountryData2
###################################################
mapParams <- mapCountryData( sPDF, nameColumnToPlot="BIODIVERSITY"
                           , addLegend=FALSE )
do.call( addMapLegend, c(mapParams, legendWidth=0.5, legendMar = 2))


###################################################
### code chunk number 10: mapGriddedData1
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
data(gridExData)
mapGriddedData(gridExData)


###################################################
### code chunk number 11: mapHalfDegreeGridToCountries1
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
mapParams <- mapHalfDegreeGridToCountries(gridExData, addLegend=FALSE)
do.call( addMapLegend, c(mapParams, legendWidth=0.5, legendMar = 2))


###################################################
### code chunk number 12: country2Region1
###################################################
#Using country2Region to calculate mean ENVHEALTH in Stern regions.
sternEnvHealth <- country2Region(inFile=countryExData
                                ,nameDataColumn="ENVHEALTH"
                                ,joinCode="ISO3"
                                ,nameJoinColumn="ISO3V10"
                                ,regionType="Stern"
                                ,FUN="mean"
                                )

print(sternEnvHealth)


###################################################
### code chunk number 13: mapByRegion1
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
mapByRegion( countryExData
           , nameDataColumn="CLIMATE"
           , joinCode="ISO3"
           , nameJoinColumn="ISO3V10"
           , regionType="Stern"
           , FUN="mean"
           )


###################################################
### code chunk number 14: finalFigure1
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")

#joining the data to a map
sPDF <- joinCountryData2Map(countryExData,
              joinCode = "ISO3",
              nameJoinColumn = "ISO3V10"
              )            
#creating a user defined colour palette
op <- palette(c('green','yellow','orange','red')) 
#find quartile breaks
cutVector <- quantile(sPDF@data[["BIODIVERSITY"]],na.rm=TRUE)
#classify the data to a factor
sPDF@data[["BIOcategories"]] <- cut(sPDF@data[["BIODIVERSITY"]]
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
              , missingCountryCol='white'
              )


###################################################
### code chunk number 15: finalFigure2
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
              , borderCol='black'
              )

## At end of plotting, reset palette to previous settings:
palette(op)  



###################################################
### code chunk number 16: bubblePlot
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")

mapBubbles( dF=getMap()
          , nameZSize="POP_EST"
          , nameZColour="GEO3major"
          , colourPalette='rainbow'
          , oceanCol='lightblue'
          , landCol='wheat'
          ) 



###################################################
### code chunk number 17: classInt_RColorBrewer
###################################################
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")

library(RColorBrewer)
library(classInt)

#getting example data and joining to a map
data("countryExData",envir=environment(),package="rworldmap")
sPDF <- joinCountryData2Map( countryExData
                           , joinCode = "ISO3"
                           , nameJoinColumn = "ISO3V10"
                           , mapResolution='coarse'
                           )

#getting class intervals using a 'jenks' classification in classInt package
classInt <- classInt::classIntervals( sPDF[["EPI"]], n=5, style="jenks")
catMethod = classInt[["brks"]]

#getting a colour scheme from the RColorBrewer package
colourPalette <- RColorBrewer::brewer.pal(5,'RdPu')

#calling mapCountryData with the parameters from classInt and RColorBrewer
mapParams <- mapCountryData( sPDF
                           , nameColumnToPlot="EPI"
                           , addLegend=FALSE
                           , catMethod = catMethod
                           , colourPalette = colourPalette )
                           
do.call( addMapLegend
       , c( mapParams
          , legendLabels="all"
          , legendWidth=0.5
          , legendIntervals="data"
          , legendMar = 2 ) )


###################################################
### code chunk number 18: multiFrame
###################################################

#uses sPDF from the previous chunk

#10 frames July 2013 started getting unableto allocate bitmap error
#op <- par(fin=c(7,9),mfcol=c(5,2),mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
#reducing to 8 frames to try to avoid memory error
op <- par(fin=c(7,9),mfcol=c(4,2),mai=c(0,0,0.2,0),xaxs="i",yaxs="i")

brewerList <- c("Greens","Greys","Oranges","OrRd"
               ,"PuBuGn","Purples","YlGn","YlGnBu","YlOrBr","YlOrRd")

#for(i in 1:10)
for(i in 1:8)
   {
    #getting a colour scheme from the RColorBrewer package
    colourPalette <- brewer.pal(7,brewerList[i])
    
    #calling mapCountryData with the parameters from RColorBrewer
    mapParams <- mapCountryData( sPDF
                               , nameColumnToPlot="CLIMATE"
                               , addLegend=FALSE
                               , colourPalette=colourPalette
                               , mapTitle=brewerList[i] )
    do.call( addMapLegend
           , c(mapParams,horizontal=FALSE,legendLabels="none",legendWidth=0.7))
   }
par(op)   


