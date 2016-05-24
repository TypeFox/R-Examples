##########################
##########################
#various pems imports
##########################
##########################

#kr 10/07/2011 

#includes 
#(functions/code below) 
##########################
#importTAB2PEMS
#importCSV2PEMS
#importOBS2PEMS
#importOB12PEMS
#

#to do
##########################
#other importers
##
#LCC GPS
#FTIR
#VBOX
#Driver Behaviour
#others as supplied
#CATI
#
##########################
#importOBS2PEMS
##
#fuel handling
#log.rate correction
#tidy code 

#comments
##########################




##########################
##########################
##import2PEMS
##########################
##########################


#kr 01/02/2012 v 0.1.0

#what it does
##########################
#import simple files to PEMS
#
#expects time.stamp and local.time to identified
#expects units to be assigned
#

#to do
##########################
#tidy code


#comments
##########################


import2PEMS <- function(file.name = file.choose(), time.stamp = NULL, local.time = NULL,
                        time.format = NULL, units = NULL, constants = NULL, history = NULL, 
                        ..., file.type = NULL, file.reader = read.delim){

    #setup
    this.call <- match.call()

    #time.format
    if(is.null(time.format))
        time.format <- "%d/%m/%Y %H:%M:%OS"


##################
#file.type code to add
##################



##################
#this needs better handling
#but can seem to document the above string
###################

    #get data

#############
#get formals of function
#use them to decide what to pass on in import
#then strip those before going on
#############


    data <- file.reader(file.name, header=TRUE)

################
#not always GMT!
################

    #sort time.stamp if there
    if(!is.null(time.stamp)){

        if(is.numeric(time.stamp)){
            data[, time.stamp] <- as.POSIXct(strptime(data[, time.stamp], format = time.format, "GMT"))
            names(data)[time.stamp] <- "time.stamp"
        }

        if(is.character(time.stamp)){
            data[, time.stamp] <- as.POSIXct(strptime(data[, time.stamp], format = time.format, "GMT"))
            names(data)[which(names(data)==time.stamp)] <- "time.stamp"
        }

    }

##################
#sort local.time if there/not there 
##################

#################
#sensible unit handler
##################

    output <- makePEMS(x = data, units = units, constants = constants, 
                       history = history, ...) 
    
    #reset history?
    if(is.null(output[["histroy"]])){
         output[["history"]] <- this.call
    } else {
         temp <- output[["history"]]
         temp[length(temp)] <- this.call
         output[["history"]] <- temp
    }

    output    


}





##########################
##########################
##import2PEMS wrappers
##########################
##########################

#kr 01/02/2012 v 0.1.0

#what it does
##########################
#import a tab delimited or clipboard file to PEMS
#import a comma delimited file to PEMS
#

importTAB2PEMS <- function(..., file.reader = read.delim) import2PEMS(..., file.reader = file.reader)
importCSV2PEMS <- function(..., file.reader = read.csv) import2PEMS(..., file.reader = file.reader)









##########################
##########################
##importOBS2PEMS
##########################
##########################

#kr 10/07/2011 v 0.3.5

#what it does
##########################
#import a standard OBS-1300 tab delimited file
##
#adds time.stamp based on reported start.time/date
#sets units
#


#to do
##########################
#tidy code
#this could be done better now
###################
#foo tidy
##v1 renamed as rename and tidy
##could do a version two and 
##put unnamed functions in sapply calls
###################
#

#comments
##########################


importOBS2PEMS <- function(file.name = file.choose(), pems = "Horiba OBS", 
          constants = NULL, history = NULL, 
          analytes = c("co", "co2", "nox", "hc"),  
          fuel = c("petrol", "diesel", "gasoline"), ...){

################################
#could do this a lot better now
################################

    #setup
    this.call <- match.call()

    #create fuel.constants
    fuel <- checkOption(fuel[1], formals(importOBS2PEMS)$fuel, 
                        "fuel", "known fuel types", 
                        fun.name = "importOBS2PEMS")
    fuel.constants <- list()
    if(fuel == "diesel")
        fuel.constants <- ref.diesel
    if(fuel == "petrol" | fuel == "gasoline")
        fuel.constants <- ref.petrol

##########################
##should this be OS not S?
##########################

    #time and date stamps
    time.stamp <- scan(file.name, nlines = 1, what = character(), quiet=TRUE)
    time.stamp <- paste(time.stamp[2], time.stamp[4], sep=" ")
    time.stamp <- as.POSIXct(strptime(time.stamp, format = "%Y/%m/%d %H:%M:%S", "GMT"))

    #read headers
    data.names <- scan(file.name, skip=1, what = character(), nlines = 1, quiet=TRUE, sep="\t")
    data.names <- gsub(" ", ".", data.names) #replace space with "." in header names
    data.names <- gsub("/", ".", data.names) #replace "/" with "." in n/s and e/w header names
    data.names <- tolower(data.names) #simplify naming
    data.names[1] <- "local.time"

    #rename analytes conc.x
    rename <- function(ans, analyte) #make analyte identifier conc.analyte so emission names are unique 
               {if(ans==analyte) ans <- paste("conc.",ans,sep="") else ans}
    for(i in 1:length(analytes)){ data.names <- sapply(data.names, rename, USE.NAMES=FALSE, analyte=analytes[i]) }

    #read units
    data.units <- scan(file.name, skip=2, what = character(), nlines = 1, quiet=TRUE, sep="\t")
    data.units <- c("Y-M-D H:M:S GMT",data.units)
    tidy<- function(ans) #strip brackets from strings
              {if(!ans=="")  
                  {if(substr(ans,1,1)=="(" & substr(ans,nchar(ans),nchar(ans))==")" ) (substr(ans,2,(nchar(ans)-1))) else ans} 
              else ans}
    data.units <- sapply(data.units, tidy, USE.NAMES=FALSE)

################
#currently gps not handled
################

    #constants
    
#this could be tidier
#move ref.chem out of this 
#then if here use else get from ref.chem?

    temp <- list(log.rate = 1000,
                 delay.co = 3.2, delay.co2 = 3.3, delay.hc = 3.9, delay.nox = 1.6, delay.afr = 1.6,
                 conc.o2 = 20.6, 
                 thc.c6 = 10, pitot.k = 2537.6, pitot.z = 0.015026, setting.gps.port = 0, k.wgec = 3.5, setting.velocity = 1, 
                 setting.720nox = 1, setting.gps = 1, setting.hc = 0, setting.option = "", setting.coco2hc = "", setting.afr = "") 
    temp[names(fuel.constants)] <- fuel.constants
    temp[names(ref.chem)] <- ref.chem
    if(is.list(constants))
        temp[names(constants)] <- constants

    extra.args <- list(...)
    temp[names(extra.args)] <- extra.args

    constants <- temp

    #read data

    data <- read.delim(file.name, header=FALSE, skip=3)

    #fix any mistmatched data and data.names sizes

    if(length(data.names)>ncol(data)){
    
        #add empties
        data[, (ncol(data)+1):length(data.names)] <- NA

    } else if(length(data.names)>ncol(data)){

        #make some missing names
        data.names <- c(data.names, paste("V", (length(data.names)+1):ncol(data), sep=""))

    }

    #add names to main data

    names(data) <- data.names

    #reset for log.rate
    #this is logger interval in ms

#tigthen this for units = seconds?
#probably not needed because it is not right in original

    data$local.time <- data$local.time * (constants$log.rate/1000)

    #pack data

    data<-cbind(time.stamp = (time.stamp + data$local.time), data)

    #gps in d.deg lat, lon
    #set hemisphere

#needs abs? reset for doing signs from scratch
#if used elsewhere

    if("latitude" %in% names(data) | "n.s" %in% names(data)){
        
        #north/south - as lower case 1 character 
        temp <- substr(tolower(as.character(data$n.s)),1,2)
        temp <- ifelse(is.na(temp), "n", temp)

        data$latitude <- ifelse(temp == "n", data$latitude, -data$latitude)
        if(!all(is.na(temp)))
            data.units[which(names(data)=="latitude")] <- "d.degLat"

    }

    if("longitude" %in% names(data) | "w.e" %in% names(data)){
        
        #east/west - as lower case 1 character 
        temp <- substr(tolower(as.character(data$w.e)),1,2)
        temp <- ifelse(is.na(temp), "e", temp)

        data$longitude <- ifelse(temp == "w", -data$longitude, data$longitude)
        if(!all(is.na(temp)))
            data.units[which(names(data)=="longitude")] <- "d.degLon"

    }

    output <- makePEMS(x = data, units = data.units, constants = constants, 
                       history = history, pems = pems, ..., silent = TRUE) 
    
    output

}









##########################
##########################
##importOB12PEMS
##########################
##########################

#kr 12/09/2013 v 0.1.0

#what it does
##########################
#import OB1 files previously created using observer
##
#adds time.stamp based on reported start.time/date
#sets units
#


#to do
##########################
#tidy code
#this could be done better now
###################
#to tidy
###################
#

#comments
##########################


importOB12PEMS <- function(file.name = file.choose(), pems = "Horiba OBS", 
          constants = NULL, history = NULL, 
          analytes = c("co", "co2", "nox", "hc"),  
          fuel = c("petrol", "diesel", "gasoline"), ...){

################################
#could do this a lot better now
################################

    #setup
    this.call <- match.call()

    #create fuel.constants
#not sure fuel matters for OB1
#already set

    fuel <- checkOption(fuel[1], formals(importOBS2PEMS)$fuel, 
                        "fuel", "known fuel types", 
                        fun.name = "importOBS2PEMS")
    fuel.constants <- list()
    if(fuel == "diesel")
        fuel.constants <- ref.diesel
    if(fuel == "petrol" | fuel == "gasoline")
        fuel.constants <- ref.petrol

    #checks
    suspect.tags <- ""

    #read first header
    input <- readLines(file.name, n = 14)
    temp <- strsplit(input[1:6], "\t")

    if(temp[[1]][2]!= "OBServer 1.0.0.")
        suspect.tags <- c(suspect.tags, "file source unrecognised")
    if(temp[[1]][3]!= "OB1 format version 1.1")
           suspect.tags <- c(suspect.tags, "suspect OB1 format")
    
    #date (yyyy/mm/dd) from line 2
    date <- as.character(temp[[2]][2])

    if(length(suspect.tags)>1)
        print(paste(suspect.tags, sep="\n"))
 
    constants <- as.list(as.numeric(temp[[6]]))
    
    #reset names 
    temp <- temp[[5]]
    temp[temp=="Log Rate"] <- "log.rate"
    temp[temp=="CO delay time"] <- "delay.co"
    temp[temp=="CO2 delay time"] <- "delay.co2"
    temp[temp=="HC delay time"] <- "delay.hc"
    temp[temp=="NOx delay time"] <- "delay.nox"
    temp[temp=="AFR delay time"] <- "delay.afr"
    temp[temp=="O2 conc"] <- "conc.o2"

    temp[temp=="CO mol mass"] <- "mm.co"
    temp[temp=="CO2 mol mass"] <- "mm.co2"
    temp[temp=="NOx mol mass"] <- "mm.nox"
    temp[temp=="C mol mass"] <- "mm.c"
    temp[temp=="H mol mass"] <- "mm.h"
    temp[temp=="O mol mass"] <- "mm.0"

    temp[temp=="PITOT k"] <- "pitot.k"
    temp[temp=="PITOT z"] <- "pitot.z"

    temp[temp=="THC/C6"] <- "thc.c6" 
    temp[temp=="Exhaust density"] <- "density.exhaust"
    temp[temp=="Fuel density"] <- "density.fuel"

    temp[temp=="COCO2HC select"] <- "setting.coco2hc"
    temp[temp=="720NOx select"] <- "setting.720nox"
    temp[temp=="K(wgec)"] <- "k.wgec" 
    temp[temp=="H/C"] <- "alpha.hc" 
    temp[temp=="Exh H/C"] <- "alpha.exhaust.hc"
    temp[temp=="GPS port"] <- "setting.gps.port"
    temp[temp=="velocity select"] <- "setting.velocity"
    temp[temp=="O/C"] <- "beta.oc"  
    temp[temp=="GPS select"] <- "setting.gps" 
    temp[temp=="HC select"] <- "setting.hc"
    temp[temp=="OPTION select"] <- "setting.option"
    temp[temp=="AFR select"] <- "setting.afr" 
    temp[temp=="Fuel select"] <- "setting.fuel"
    temp[temp=="Vehicle type"] <- "vehicle.type"

    names(constants) <- temp

    test.fun <- function(x)
                     if(is.na(names(constants)[x])  || names(constants)[x]=="O") 
                         FALSE else TRUE

    constants <- constants[sapply(1:length(constants), test.fun)]

    history <- list(input, this.call)

    units <- read.delim(file.name, skip=14, nrows=1, header=TRUE, stringsAsFactors=FALSE)

    names(units) <- tolower(names(units))
    names(units) <- gsub("real.time.fuel.consumption.by.", "rtfc.", names(units))
    names(units) <- gsub("fuel.consumption.by.", "fc.", names(units))
    names(units) <- gsub("afr.by.", "afr.", names(units))
    names(units) <- gsub("option[.]", "option", names(units))
    names(units) <- gsub(".sensor", "", names(units))
    names(units) <- gsub("pitot.output..retro.calc.", "pitot.output", names(units))

    names(units)[names(units)=="time.1"] <- "time.stamp"
    names(units)[names(units)=="time"] <- "local.time"
    for(i in analytes){
        names(units)[names(units)==i]<- paste("conc.", i, sep="")
        names(units)[names(units)==paste(i, ".1", sep="")]<- paste("em.", i, sep="")
    }
    
    temp <- as.character(units)
    temp <- gsub("[(]", "", temp)
    temp <- gsub("[)]", "", temp)
    temp[temp=="hh:mm:ss"] <- "Y-M-D H:M:S GMT"
    temp[temp=="NA"] <- ""
    temp[temp=="g/sec"] <- "g/s"
    units[] <- temp 

    data <- read.delim(file.name, skip=16, header=FALSE, na.strings = c("NA", "Inf (No Carbon)"))
    names(data) <- names(units)
    data$time.stamp <- paste(date, data$time.stamp, sep=" ")

    #earlier versions of OB1 files did not set lat lon sign
    #later versions did 
    #so cannot assume systematic handling
    #so reset all to positive and redo signs based on n.s and w.e
    
    if("latitude" %in% names(data) & "n.s" %in% names(data)){
        
        #north/south - as lower case 1 character 
        temp <- substr(tolower(as.character(data$n.s)),1,2)
        temp <- ifelse(is.na(temp), "n", temp)

        data$latitude <- abs(data$latitude)
        data$latitude <- ifelse(temp == "n", data$latitude, -data$latitude)

        temp[names(units)=="latitude"] <- "d.degLat"
    }

    if("longitude" %in% names(data) & "w.e" %in% names(data)){
        
        #east/west - as lower case 1 character 
        temp <- substr(tolower(as.character(data$w.e)),1,2)
        temp <- ifelse(is.na(temp), "e", temp)

        data$longitude <- abs(data$longitude)
        data$longitude <- ifelse(temp == "w", -data$longitude, data$longitude)
        
        temp[names(units)=="longitude"] <- "d.degLon"
    }

    output <- makePEMS(x = data, units = units, constants = constants, 
                       history = history, pems = pems, ..., silent=TRUE) 
    
    output
    }
