### R code from vignette source 'browseMovebank.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("move")


###################################################
### code chunk number 2: createCURLHandle
###################################################
#curl <- movebankLogin(username="user", password="password")


###################################################
### code chunk number 3: seatchMovebankStudies
###################################################
#searchMovebankStudies(x="oose", login=curl)


###################################################
### code chunk number 4: getMovebankID
###################################################
#getMovebankID("BCI Ocelot",login=curl)
#> 123413


###################################################
### code chunk number 5: getMovebankStudy
###################################################
#getMovebankStudy(study="BCI Ocelot",login=curl)


###################################################
### code chunk number 6: getMovebankSensors
###################################################
#getMovebankSensors("BCI Ocelot",login=curl)


###################################################
### code chunk number 7: getMovebankSensorsAttributes
###################################################
#getMovebankSensorsAttributes("BCI Ocelot",login=curl)


###################################################
### code chunk number 8: getMovebankAnimals
###################################################
#getMovebankAnimals(study="BCI Ocelot",login=curl)


###################################################
### code chunk number 9: getMovebankData
###################################################
###Download timestamp and coordinates
#bobby <- getMovebankData(study="BCI Ocelot", animalName="Bobby", login=curl)

###Create a move object
#bobby <- getMovebankData(study="BCI Ocelot", animalName="Bobby", 
#                                               login=curl, moveObject=TRUE)


