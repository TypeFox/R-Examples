################################
##
## Class: Dataclass
##
################################


## Access methods
if(!isGeneric("filename")) 
   setGeneric("filename", function(object) standardGeneric("filename"))
if(!isGeneric("Data")) 
   setGeneric("Data", function(object) standardGeneric("Data"))
if(!isGeneric("runs")) 
   setGeneric("runs", function(object) standardGeneric("runs"))
if(!isGeneric("obsDim"))
   setGeneric("obsDim", 
               function(object) standardGeneric("obsDim")) ### new v.1.8
if(!isGeneric("name")) 
   setGeneric("name", function(object) standardGeneric("name"))### new v.1.8
if(!isGeneric("getVersion")) 
  setGeneric("getVersion", 
              function(object) standardGeneric("getVersion"))### new v.1.8

# Replacement methods
if(!isGeneric("filename<-")) 
   setGeneric("filename<-", 
               function(object, value) standardGeneric("filename<-"))
if(!isGeneric("Data<-")) 
   setGeneric("Data<-", function(object, value) standardGeneric("Data<-"))
if(!isGeneric("name<-")) 
   setGeneric("name<-", 
                function(object, value) standardGeneric("name<-")) ### new v.1.8

### not be set for Dataclass itself without generating a new object:

if(!isGeneric("runs<-")) 
   setGeneric("runs<-", function(object, value) standardGeneric("runs<-"))
if(!isGeneric("samplesize<-")) 
   setGeneric("samplesize<-", 
               function(object, value) standardGeneric("samplesize<-"))

### new 2.0
if(!isGeneric("Dataclass")) 
    setGeneric("Dataclass", 
              function(Data,...) standardGeneric("Dataclass"))

rbind <- function(x, ...) base::rbind(x,...)

if(!isGeneric("rbind")) 
    setGeneric("rbind", 
              function(x,...) standardGeneric("rbind"))

if(!isGeneric("obsdimnames")) 
    setGeneric("obsdimnames", function(object) standardGeneric("obsdimnames"))
if(!isGeneric("obsdimnames<-")) 
    setGeneric("obsdimnames<-", function(object, value) standardGeneric("obsdimnames<-"))
if(!isGeneric("runnames")) 
    setGeneric("runnames", function(object) standardGeneric("runnames"))
if(!isGeneric("runnames<-")) 
    setGeneric("runnames<-", function(object, value) standardGeneric("runnames<-"))


# general methods

## moved to distr version 1.9
#if(!isGeneric("isOldVersion")) 
#   setGeneric("isOldVersion", function(object) standardGeneric("isOldVersion"))
#
#if(!isGeneric("conv2NewVersion")) 
#   setGeneric("conv2NewVersion", 
#               function(object) standardGeneric("conv2NewVersion"))

if(!isGeneric("savedata")) 
    setGeneric("savedata", function(object,...) standardGeneric("savedata"))

################################
##
## Class: Simulation
##
################################

### changed from version 1.8 on:
## ith observation in ith line of datamatrix/array
## jth item/dimension of each observation in jth column of datamatrix/array
## kth run/time of each observation in kth slide of datamatrix/array

## ++old
## +ith run in ith line of datamatrix
## +jth samples of each run in jth column of datamatrix



## Access Methods
if(!isGeneric("seed")) 
   setGeneric("seed", function(object) standardGeneric("seed"))

## Replace Methoden
if(!isGeneric("seed<-")) 
   setGeneric("seed<-", function(object, value) standardGeneric("seed<-"))
if(!isGeneric("distribution<-")) 
   setGeneric("distribution<-", 
                function(object, value) standardGeneric("distribution<-"))

## general methods

## Simulation method
if(!isGeneric("simulate")) 
   setGeneric("simulate",
               function(object, nsim=-1, seed=-1, ...)
                         standardGeneric("simulate"))

################################
##
## Class: Contsimulation
##
################################

## Access methods
if(!isGeneric("ind")) 
   setGeneric("ind", function(object) standardGeneric("ind"))
if(!isGeneric("Data.id")) 
   setGeneric("Data.id", function(object) standardGeneric("Data.id"))
if(!isGeneric("Data.c")) 
    setGeneric("Data.c", function(object) standardGeneric("Data.c"))
if(!isGeneric("rate")) 
    setGeneric("rate", function(object) standardGeneric("rate"))
if(!isGeneric("distribution.c")) 
   setGeneric("distribution.c", 
               function(object) standardGeneric("distribution.c"))
if(!isGeneric("distribution.id")) 
   setGeneric("distribution.id", 
               function(object) standardGeneric("distribution.id"))
if(!isGeneric("seed")) 
    setGeneric("seed", function(object) standardGeneric("seed"))


## Replace Methoden

##if(!isGeneric("rate<-")) 
#    setGeneric("rate<-", function(object, value) standardGeneric("rate<-"))
#if(!isGeneric("distribution.c<-")) 
    setGeneric("distribution.c<-", 
                function(object, value) standardGeneric("distribution.c<-"))
#if(!isGeneric("distribution.id<-")) 
    setGeneric("distribution.id<-", 
                function(object, value) standardGeneric("distribution.id<-"))
#if(!isGeneric("seed<-")) 
    setGeneric("seed<-", function(object, value) standardGeneric("seed<-"))
