## ----install, eval=FALSE-------------------------------------------------
#  library(devtools)
#  install_github("hackout2/repijson")

## ----load----------------------------------------------------------------
library("repijson")

## ---- echo=FALSE, fig.width=7, fig.height=5, fig.path="figs/"------------
#here deliberately not echoed just to show the diagram
epijsonObjectVis(textSize = 5) #textSize 4 is default

## ---- echo=FALSE, fig.width=5, fig.height=5, fig.path="figs/"------------
#here deliberately not echoed just to show the diagram
epijsonObjectVis( attribMeta = 'ejAttribute',
                  attribRecord = 'ejAttribute',
                  attribEvent = 'ejAttribute',
                  labelObject = 'ejObject',
                  labelMeta = 'ejMetadata',
                  labelRecord = 'ejRecord',
                  labelEvent = 'ejEvent',
                  textSize = 5) #textSize 4 is default

## ------------------------------------------------------------------------
toyll[1:3,1:8]

## ------------------------------------------------------------------------
#converting dates to date format
toyll$date.of.admission <- as.POSIXct(toyll$date.of.admission)
toyll$date.of.discharge <- as.POSIXct(toyll$date.of.discharge)
#create ejObject
ejOb <- as.ejObject(toyll[1:3,],
                  recordAttributes=c("name","gender"),
                  eventDefinitions=list(
                      define_ejEvent(name="admission", date="date.of.admission",attributes="hospital"),
                      define_ejEvent(name="discharge", date="date.of.discharge")
                  ))
#display ejObject
ejOb

## ------------------------------------------------------------------------
library(OutbreakTools)
library(sp)
library(HistData)

## ------------------------------------------------------------------------
data(Snow.deaths)

## ------------------------------------------------------------------------
simulated <- Snow.deaths
simulated$gender <- c("male","female")[(runif(nrow(simulated))>0.5) +1]
simulated$date <- as.POSIXct("1854-04-05") + rnorm(nrow(simulated), 10) * 86400
simulated$pump <- ceiling(runif(nrow(simulated)) * 5)

exampledata1<-head(simulated)
exampledata1

## ------------------------------------------------------------------------
exampledata2<- data.frame(id=c(1,2,3,4,5),
                 name=c("Tom","Andy","Ellie","Ana","Tibo"),
                 dob=c("1981-01-12","1980-11-11","1982-02-10","1981-12-09","1983-03-08"),
                 gender=c("male","male","female","female","male"),
                 date.of.onset=c("2014-12-28","2014-12-29","2015-01-03","2015-01-08","2015-01-04"),
                 date.of.admission=c(NA,"2015-01-05","2015-01-12",NA,"2015-01-14"),
                 date.of.discharge=c(NA,NA,"2015-01-17",NA,"2015-01-17"),
                 hospital=c(NA,"St Marys","Whittington",NA,"Whittington"),
                 fever=c("yes","yes","no","no","yes"),
                 sleepy=c("no","yes","yes","no","yes"),
                 contact1.id=c("B","A","5",NA,"3D"),
                 contact1.date=c("2014-12-26","2014-12-26","2014-12-28",NA,"2014-12-28"),
                 contact2.id=c("3D","3D","5",NA,"A"),
                 contact2.date=c("2014-12-25","2014-12-26","2015-01-14",NA,"2014-12-25"),
                 contact3.id=c("B",NA,NA,NA,NA),
                 contact3.date=c("2015-01-08",NA,NA,NA,NA)
                 )
            
exampledata2

## ------------------------------------------------------------------------
eg1 <- as.ejObject(exampledata1,	
    recordAttributes = "gender",	
    eventDefinitions = list(define_ejEvent(date="date",	name="Death", location=list(x="x", y="y", proj4string=""), attributes="pump")),
 		metadata=list())		       
eg1

## ------------------------------------------------------------------------
exampledata2$date.of.onset <- as.POSIXct(exampledata2$date.of.onset)
exampledata2$date.of.admission <- as.POSIXct(exampledata2$date.of.admission)
exampledata2$date.of.discharge <- as.POSIXct(exampledata2$date.of.discharge)
exampledata2$contact1.date <- as.POSIXct(exampledata2$contact1.date)
exampledata2$contact2.date <- as.POSIXct(exampledata2$contact2.date)
exampledata2$contact3.date <- as.POSIXct(exampledata2$contact3.date)

## ------------------------------------------------------------------------
eg2 <- as.ejObject(exampledata2, recordAttributes = c("id","name","dob","gender"),
     eventDefinitions = list(define_ejEvent(name="Date Of Onset", date="date.of.onset", 
                                            attributes=list()),
                             define_ejEvent(name="Hospital admission", date="date.of.admission", 
											attributes=list("hospital", "fever", "sleepy")),
							 define_ejEvent(name="Hospital discharge", date="date.of.discharge"),
							 define_ejEvent(name="Contact1", date="contact1.date", attributes=list("contact1.id")),
							 define_ejEvent(name="Contact2", date="contact2.date", attributes=list("contact2.id")),
							 define_ejEvent(name="Contact3", date="contact3.date", attributes=list("contact3.id"))
						),
 		metadata=list())
eg2

## ------------------------------------------------------------------------
as.data.frame(eg1)
as.data.frame(eg2)

## ------------------------------------------------------------------------
data(ToyOutbreak) 

## ------------------------------------------------------------------------
eg3 <- as.ejObject(ToyOutbreak)

## ---- fig.path="figs/"---------------------------------------------------
sp_eg1 <- as.SpatialPointsDataFrame.ejObject(eg1)
plot(sp_eg1,pch=20,col="green")
text(10,17,"Example from Snow Deaths data")

