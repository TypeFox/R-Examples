library(micromap)

# Code for Section 1 Introduction
### intial example
data("edPov")
head(edPov)


data("USstates")
head(USstates@data)


statePolys <- create_map_table(USstates, 'ST')
head(statePolys)


### draft micromap plot figure 1
mmplot(stat.data=edPov,
       map.data=statePolys,
       panel.types=c('labels', 'dot', 'dot','map'),
       panel.data=list('state','pov','ed', NA),
       ord.by='pov',   
       grouping=5, median.row=T,
       map.link=c('StateAb','ID'))


### publication figure 2
mmplot(stat.data=edPov,  map.data=statePolys ,
       panel.types=c('labels', 'dot', 'dot','map'),
       panel.data=list('state','pov','ed', NA),
       ord.by='pov',  
       grouping=5, 
       median.row=T,
       map.link=c('StateAb','ID'),
       
       plot.height=9,							
       colors=c('red','orange','green','blue','purple'), 
       map.color2='lightgray',
       
       
       panel.att=list(list(1, header='States', panel.width=.8, align='left', text.size=.9),
                      list(2, header='Percent Living Below \n Poverty Level',
                           graph.bgcolor='lightgray', point.size=1.5,
                           xaxis.ticks=list(10,15,20), xaxis.labels=list(10,15,20),
                           xaxis.title='Percent'),
                      list(3, header='Percent Adults With\n4+ Years of College',
                           graph.bgcolor='lightgray', point.size=1.5,
                           xaxis.ticks=list(0,20,30,40), xaxis.labels=list(0,20,30,40),
                           xaxis.title='Percent'),
                      list(4, header='Light Gray Means\nHighlighted Above',  
                           inactive.border.color=gray(.7), inactive.border.size=2,	
                           panel.width=.8)))



# Code for Section 2 Quick Plotting tips
### publication figure 3
mmplot (stat.data=edPov, map.data=statePolys,
        panel.types=c('dot_legend', 'labels', 'dot', 'dot', 'map'),
        panel.data=list(NA, 'state', 'pov', 'ed', NA),
        map.link=c('StateAb','ID'),
        ord.by='pov', 
        grouping=5, 
        median.row=T,
                
        plot.height=9, 
        
        colors=c('red','orange','green','blue','purple'),
        map.color2='lightgray', 
        
        panel.att=list(list(1, point.type=20, point.size=2, point.border=TRUE,
                            graph.border.color='white',
                            xaxis.text.display=F, xaxis.line.display=F,
                            graph.grid.major=F),
                       
                       list(2, header='States', panel.width=.8, 
                            align='left', text.size=.9),
                       
                       list(3, header='Percent Living Below\nPoverty Level',
                            graph.bgcolor='lightgray', point.size=1.5,
                            xaxis.ticks=list(10,15,20), xaxis.labels=list(10,15,20),
                            xaxis.title='Percent'),
                       
                       list(4, header='Percent Adults With\n4+ Years of College',
                            graph.bgcolor='lightgray', point.size=1.5,
                            xaxis.ticks=list(20,30,40), xaxis.labels=list(20,30,40), 
                            xaxis.title='Percent', left.margin=-.8, right.margin=0),
                       
                       list(5, header='Light Gray Means\nHighlighted Above', 
                            inactive.border.color=gray(.7), inactive.border.size=2, 
                            panel.width=.8)))




### publication figure 4
myPlot <- mmplot(stat.data=edPov, map.data=statePolys,
                 panel.types=c('map', 'dot_legend',  'labels', 'dot', 'dot'),
                 panel.data=list(NA, NA, 'state', 'pov', 'ed'),
                 map.link=c('StateAb','ID'),
                 ord.by='pov', 
                 grouping=5, 
                 median.row=T,
                                  
                 plot.height=9,
                 
                 colors=c('red','orange','green','blue','purple'),
                 map.color2='lightgray', 
                 
                 panel.att=list(list(2, point.type=20, point.size=2, point.border=TRUE,
                                     graph.border.color='white',
                                     xaxis.text.display=F, xaxis.line.display=F,
                                     graph.grid.major=F),
                                
                                list(3, header='States', panel.width=.8, 
                                     align='left', text.size=.9),
                                
                                list(4, header='Percent Living Below\nPoverty Level',
                                     graph.bgcolor='lightgray', point.size=1.5,
                                     xaxis.ticks=list(10,15,20), xaxis.labels=list(10,15,20),
                                     xaxis.title='Percent'),
                                
                                list(5, header='Percent Adults With\n4+ Years of College',
                                     graph.bgcolor='lightgray', point.size=1.5,
                                     xaxis.ticks=list(20,30,40), xaxis.labels=list(20,30,40), 
                                     xaxis.title='Percent'),
                                
                                list(1, header='Light Gray Means\nHighlighted Above', 
                                     inactive.border.color=gray(.7), inactive.border.size=2, 
                                     panel.width=.8)))

print(myPlot, name='myExhibit.tiff', res=300)



# Code for Section 3 Preparing data for use with the library
# The download.file command may or may not work depending on firewall settings.
# If it fails, then download the file from the ftp directory to a local
# directory manually.
download.file("ftp://ftp.epa.gov/wed/ecoregions/tx/tx_eco_l3.zip",
              "C:/temp/tx_eco_l3.zip")
library(rgdal)
library(utils)
unzip("C:/temp/tx_eco_l3.zip", exdir="C:/temp")
txeco <- readOGR("C:/temp","tx_eco_l3")
plot(txeco)

# Create an ID column in your spatial dataframe for the 
# create_map_table function
txeco$ID <- txeco$US_L3CODE
# Load the maptools library
library(maptools)
txeco1<-thinnedSpatialPoly(txeco,tolerance=7000, minarea=0.001, topologyPreserve = F, avoidGEOS = T)
plot(txeco1)
# Alternatively, load the rgeos library
library(rgeos)
txeco2 <- gSimplify(txeco, 7000, topologyPreserve=T)
class(txeco2)
txeco2 <- SpatialPolygonsDataFrame(txeco2, txeco@data)
class(txeco2)
plot(txeco2)

# Simplifying very large spatial data sets
# The download.file command may or may not work depending on firewall settings.
# If it fails, then download the file from the ftp directory to a local
# directory manually.
download.file("ftp://ftp.epa.gov/wed/ecoregions/us/Eco_Level_III_US.zip",
              "C:/temp/Eco_Level_III_US.zip")
              
unzip("C:/temp/Eco_Level_III_US.zip", exdir="C:/temp")
eco3 <- readOGR("C:/temp", "us_eco_l3")
plot(eco3)              

eco3_thin1 <- thinnedSpatialPoly(eco3, tolerance=50000, topologyPreserve=TRUE, avoidGEOS=FALSE)
              
eco3_thin2 <- thinnedSpatialPoly(eco3, tolerance=50000, minarea=100,avoidGEOS= TRUE)
              
eco3_thin3 <- gSimplify(eco3, tol=50000, topologyPreserve=TRUE)
class(eco3_thin3) #note gSimplify returns a spatial polygon object

# convert spatial polygon object to spatial polygon data frame
df <- eco3@data
eco3 <- SpatialPolygonsDataFrame(eco3_thin3, df)
              
# Code for section 5 Creating a new panel type
### new graph type instructions
data("lungMort")
myStats <- lungMort
head(myStats)


myStats <- subset(myStats, !StateAb=='DC')


myNewStats <- create_DF_rank(myStats, ord.by="Rate_00", group=5)
head(myNewStats)
# step 1.
### ggplot2 code:
ggplot(myNewStats) + 
               geom_segment(aes(x=Rate_95, y=-pGrpOrd,
               xend=Rate_00, yend=-pGrpOrd, colour=factor(color)),
               arrow=arrow(length=unit(0.1,"cm"))) +
               facet_grid(pGrp~., scales="free_y") +
               scale_colour_manual(values=c('red','orange','green','blue',
               'purple'), guide="none")



# step 2.
myAtts <- sample_att()
myNumber <- 1

myAtts$colors  <- c('red','orange','green','blue','purple')
myAtts[[myNumber]]$panel.data <- c('Rate_95','Rate_00')


myColors <- myAtts$colors # pulls color out of the plot level section of the 'myAtts'attributes list
myColumns <- myAtts[[myNumber]]$panel.data 	# looks in the panel level section numbered 'myNumber' of the 'myAtts' attributes list 

myNewStats$data1 <- myNewStats[, myColumns[1] ]
myNewStats$data2 <- myNewStats[, myColumns[2] ]


## build graph
myPanel  <- ggplot(myNewStats) +
  geom_segment(aes(x=data1, y=-pGrpOrd,
                   xend= data2, yend=-pGrpOrd, colour=factor(color)),
               arrow=arrow(length=unit(0.1,"cm"))) +
  facet_grid(pGrp~.) +
  scale_colour_manual(values=myColors, 
                      guide="none")

myPanel
class(myPanel)

# step 3.
## add in plot attributes
assimilatePlot(myPanel, myNumber, myAtts) 


## build function
arrow_plot_build <- function(myPanel, myNumber, myNewStats, myAtts){
  myColors <- myAtts$colors 			
  myColumns <- myAtts[[myNumber]]$panel.data 	
  
  myNewStats$data1 <- myNewStats[, myColumns[1] ]
  myNewStats$data2 <- myNewStats[, myColumns[2] ]
  
  myPanel  <- ggplot(myNewStats) +
    geom_segment(aes(x=data1, y=-pGrpOrd,
                     xend= data2, yend=-pGrpOrd,
                     colour=factor(color)),
                 arrow=arrow(length=unit(0.1,"cm")))  +
    facet_grid(pGrp~.) +
    scale_colour_manual(values=myColors, guide="none")
  
  myPanel <- assimilatePlot(myPanel, myNumber, myAtts) 
  
  myPanel
}


## add median row code
arrow_plot_build <- function(myPanel, myNumber, myNewStats, myAtts){
  myColors <- myAtts$colors
  myColumns <- myAtts[[myNumber]]$panel.data
  myNewStats$data1 <- myNewStats[, myColumns[1]]
  myNewStats$data2 <- myNewStats[, myColumns[2]]
  myNewStats <- alterForMedian(myNewStats, myAtts)
  myPanel <- ggplot(myNewStats) +29+ geom_segment(aes(x=data1, y=-pGrpOrd,
                     xend= data2, yend=-pGrpOrd,
                     colour=factor(color)),
                     arrow=arrow(length=unit(0.1,'cm'))) +
    facet_grid(pGrp~., space="free", scales='free_y') +
    scale_colour_manual(values=myColors, guide='none')
  myPanel <- assimilatePlot(myPanel, myNumber, myAtts)
  }
myPanel

# step 4.
## build specialized attributes list
print(myAtts) #see full list of attributes available for alteration / specification by a user
myPanelAtts <- standard_att()
myPanelAtts <- append(myPanelAtts, 
                      list(line.width=1, tip.length=1))
myPanelAtts


arrow_plot_att <- function(){
  myPanelAtts <- standard_att()
  myPanelAtts <- append(myPanelAtts, 
                        list(line.width=1, tip.length=1))
}
myPanelAtts


## implement specialized attributes
arrow_plot_build <- function(myPanel, myNumber, myNewStats, myAtts){
  myColors <- myAtts$colors 			
  myColumns <- myAtts[[myNumber]]$panel.data 	
  myLineWidth <- myAtts[[myNumber]]$line.width	# Again, note that these are stored in the panel level section of the 
  myTipLength <- myAtts[[myNumber]]$tip.length	#    attributes object
  
  
  myNewStats$data1 <- myNewStats[, myColumns[1] ]
  myNewStats$data2 <- myNewStats[, myColumns[2] ]
  
  myNewStats <- alterForMedian(myNewStats, myAtts)
  
  myPanel  <- ggplot(myNewStats) +
    geom_segment(aes(x=data1, y=-pGrpOrd,
                     xend= data2, yend=-pGrpOrd, 
                     colour=factor(color)),
                 arrow=arrow(length=unit(0.1*myTipLength,"cm")),	# Here, you'll notice the '1' default above # is specifying length in tenths of a cm
                 size=myLineWidth)  +
    facet_grid(pGrp~., space='free', scales='free_y') +
    scale_colour_manual(values=myColors, guide='none')
  
  myPanel <- assimilatePlot(myPanel, myNumber, myAtts) 
}
myPanel


## figure 7 Basic micromap plot with arrow panel
data("USstates")
statePolys <- create_map_table(USstates, 'ST')
mmplot(stat.data=myStats,
       map.data=statePolys,
       panel.types=c('map','labels', 'arrow_plot'),
       panel.data=list(NA,'State', list('Rate_95','Rate_00')), 
       ord.by='Rate_00', 
       grouping=5,
       map.link=c('StateAb','ID'),
       panel.att=list(list(3, line.width=1.25, tip.length=1.5)))


## figure 8 Micromap plot with arrow panel
myStats <- lungMort

mmplot(stat.data=myStats,
       map.data=statePolys,
       panel.types=c('map', 'dot_legend', 'labels', 'dot_cl', 'arrow_plot'),
       panel.data=list(NA,
                         'points',
                         'State',
                         list('Rate_00','Lower_00','Upper_00'),
                         list('Rate_95','Rate_00')),
       ord.by='Rate_00', grouping=5,
       median.row=T,
       map.link=c('StateAb','ID'),
       plot.height=10,
       colors=c('red','orange','green','blue','purple'),
       panel.att=list(list(1, header='Light Gray Means\n Highlighted Above',
                             map.all=TRUE,
                             fill.regions='two ended',
                             inactive.fill='lightgray',
                             inactive.border.color=gray(.7),
                             inactive.border.size=2,
                             panel.width=1),
                        list(2, point.type=20,
                               point.border=TRUE),
                        list(3, header='U.S. \nStates ',
                               panel.width=.8,
                               align='left', text.size=.9),
                        list(4, header='State 2000\n Rate and 95% CI',
                               graph.bgcolor='lightgray',
                               xaxis.ticks=c(20,30,40,50),
                               xaxis.labels=c(20,30,40,50),
                               xaxis.title='Deaths per 100,000'),
                        list(5, header='State Rate Change\n 1995-99 to 2000-04',
                               line.width=1.25, tip.length=1.5,
                               graph.bgcolor='lightgray',
                               xaxis.ticks=c(20,30,40,50),
                               xaxis.labels=c(20,30,40,50),
                               xaxis.title='Deaths per 100,000'))) 


# Code for section 6 Group-Categorized Micromaps
### mmgroupedPlot
data('vegCov')
data('WSA3')
print(vegCov)
print(WSA3@data)

wsa.polys <- create_map_table(WSA3)
head(wsa.polys)

# create a national polygon area based on the 9 NARS reporting regions
national.polys <-subset(wsa.polys, hole==0 & plug==0)
national.polys <- transform(national.polys, ID='National', region=4, poly=region*1000 + poly)
head(national.polys)

wsa.polys <- rbind(wsa.polys, national.polys)


### figure 9 Basic group-categorized micromap plot
mmgroupedplot(stat.data=vegCov,
              map.data=wsa.polys,
              panel.types=c('map', 'labels', 'bar_cl', 'bar_cl'),
              panel.data=list(NA,'Category',
                              list('Estimate.P','LCB95Pct.P','UCB95Pct.P'),
                              list('Estimate.U','LCB95Pct.U','UCB95Pct.U')),
              grp.by='Subpopulation',
              cat='Category',
              map.link=c('Subpopulation', 'ID'))



### figure 10 Group-categorized micromap plot
mmgroupedplot(stat.data= vegCov,
              map.data= wsa.polys,
              panel.types=c('map', 'labels', 'bar_cl', 'bar_cl'),
              panel.data=list(NA,'Category',
                              list('Estimate.P','LCB95Pct.P','UCB95Pct.P'),
                              list('Estimate.U','LCB95Pct.U','UCB95Pct.U')),
              grp.by='Subpopulation',
              cat='Category',
              colors=c('red3','green3','lightblue'),
              map.link=c('Subpopulation', 'ID'),
              map.color='orange3',
              plot.grp.spacing=2,
              plot.width=7,
              plot.height=4,
              
              panel.att=list(list(1, header='Region', header.size=1.5, 
                                  panel.width=.75), 
                             list(2, header='Category', 
                                  header.size=1.5, 
                                  panel.width=1.7),
                             list(3, header='Percent', header.size=1.5, 
                                  graph.bgcolor='lightgray',
                                  xaxis.title='percent',
                                  xaxis.ticks=c(0,20,40,60),
                                  xaxis.labels=c(0,20,40,60)),
                             list(4, header='Unit', header.size=1.5, 
                                  graph.bgcolor='lightgray',
                                  xaxis.title='thousands',
                                  xaxis.ticks=c(0, 200000,350000,550000),
                                  xaxis.labels=c(0, 200,350,550))))
