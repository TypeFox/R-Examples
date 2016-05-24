##########################
##########################
#various pems imports
##########################
##########################


#additional (2)

#work in progress


#############################
#############################
##importParSYNC2PEMS
#############################
#############################

#quick importer for parSYNC data
#started 15-03-2015
#kr v0.0.1 (15-03-2015)

#added importer for CAGE data
#kr v0.0.1 (12-08-2015)

#two functions
###################
#importParSYNC2PEMS
#importCAGE2PEMS
#

#currently very crude because this is in-development data type


###################
#to do
###################
#everything
#


###################
#notes
###################
#

###################
#importParSYNC2PEMS
###################


importParSYNC2PEMS <- function(file.name = file.choose(), reset.signals = TRUE, 
         history = NULL, constants = NULL, pm.analyzer = "parSYNC",
         ... ){
   
    #setup
    this.call <- match.call()

#############
#time zone might need better handling
#############

##############
#not sure I need this
##############
    fun.name <- "importParSYNC2PEMS"
    extra.args <- list(...)
    to.lower <- if("to.lower" %in% names(extra.args))
                    extra.args$to.lower else TRUE
    extra.args <- extra.args[!names(extra.args) %in% "to.lower"]

    #set up import
    extra.args <- listUpdate(list(header=TRUE), extra.args)
    extra.args$file <- file.name
    ans <- do.call(read.csv, extra.args)

    #reset time stamps
    temp <- ans[c("Timestamp", "Date", "Time")]
    ans <- ans[names(ans)[!names(ans) %in% c("Timestamp", "Date", "Time")]]

#####################
#error handling?
#####################
    #stop is not valid
    ##if(is.null(ans))
        ##stop(paste("\t In ", fun.name,"(...) selected file type not recognised", sep=""),
             ##paste("\n\t [suggest setting file.type in call if known]", sep=""), 
             ##call. = FALSE, domain = NA)

    temp$time.stamp <- paste(temp$Date, temp$Time, sep=" ")


#########################
#need to check this out
#re GMT handling
#########################

    #old lines
    #first version, remove GMT default
    #    temp$time.stamp <- as.POSIXct(strptime(temp$time.stamp, format = "%m/%d/%Y %H:%M:%OS"), tz="GMT")
    #    temp$time.stamp <- as.POSIXct(strptime(temp$time.stamp, format = "%m/%d/%Y %H:%M:%OS"))

    temp$time.stamp <- if(length(grep("PM|AM", temp$time.stamp))>0){
                           as.POSIXct(strptime(temp$time.stamp, format = "%m/%d/%Y %I:%M:%S %p"))
                       } else {
                           as.POSIXct(strptime(temp$time.stamp, format = "%m/%d/%Y %H:%M:%OS"))
                       }    
    temp$local.time <- as.numeric(temp$time.stamp - temp$time.stamp[1])

    #check for dates end/clocking issue
    if(any(diff(temp$local.time)<0))
            warning("possible clocking issue with time stamp")

    temp$parsync.timestamp <- temp$Timestamp

    temp <- temp[c("time.stamp", "local.time", "parsync.timestamp")]
    ans <- cbind(temp, ans)

    #create units
    units <- rep("", length(names(ans)))
    units[grep("..V.", names(ans))] <- "V"
    units[grep("..deg.C.", names(ans))] <- "degC"
    #old line
    #   units[1] <- "Y-M-D H:M:S GMT"
    units[1] <- "Y-M-D H:M:S"
    units[2] <- "s"

    #tidy names
    names(ans) <- gsub("..V.", "", names(ans))
    names(ans) <- gsub("..deg.C.", "", names(ans))
    names(ans) <- gsub("Bag.", "Bag", names(ans))

###########################
#special handling
###########################

    #special handling for signals
    #reverse signals
    if(is.logical(reset.signals) && reset.signals){
        reset.signals <- c("Opacity", "Ionization")
    }
    if(is.character(reset.signals)){
        for(i in 1:length(reset.signals))
            if(reset.signals[i] %in% names(ans))
                ans[,reset.signals[i]] <- -ans[,reset.signals[i]]
    }

############################
#old
#    if(is.logical(reset.opacity))
#        ans$Opacity <- -ans$Opacity
##############

    if(to.lower)
        names(ans) <- tolower(names(ans))

    #make pems
    output <- makePEMS(x = ans, units = units, constants = constants, history = history, pm.analyzer=pm.analyzer, ...)
 
    class(output) <- "not.pems"
    output$history[length(output$history)] <- this.call 
    class(output) <- "pems"
   
    #return output
    return(output)



}




###################
#importCAGE2PEMS
###################


importCAGE2PEMS <- function(..., calibrator = "CAGE"){

    #setup
    this.call <- match.call()

    #currently uses parSYNC import function
    #this might change 
    ans <- importParSYNC2PEMS(...)


    #might turn off some of the 'to neg' options later
    #reset.signal = FALSE in above
    #why if not likely to be there anyway..?

    #tidy names re move parsyn to cage
    names(ans) <- gsub("parsync", "cage", names(ans))
    if("x.rh" %in% names(ans))
        names(ans)[names(ans)=="x.rh"] <- "rh"

    #tidy units
    #need tolower in parsync importer already changed these
    units(ans)[grep("..v.", tolower(names(ans)))] <- "V"
    units(ans)[grep("..deg.c.", tolower(names(ans)))] <- "degC"
    units(ans)[names(units(ans)) == "rh"] <- "%"
    units(ans)[grep("..lpm.", tolower(names(ans)))] <- "L/min"

    #tidy names re units
    names(ans) <- gsub("..V.", "", names(ans))
    names(ans) <- gsub("..v.", "", names(ans))
    names(ans) <- gsub("..deg.C.", "", names(ans))
    names(ans) <- gsub("..deg.c.", "", names(ans))
    names(ans) <- gsub("Bag.", "Bag", names(ans))
    names(ans) <- gsub("bag.", "bag", names(ans))
    names(ans) <- gsub("..lpm.", "", names(ans))

    #update tags re move parsyn to cage
    temp <- class(ans)
    class(ans) <- "break"
    ans <- ans[!names(ans) %in% "pm.analyzer"]
    ans$calibrator <- calibrator
    ans$history <- this.call 
    class(ans) <- temp

    #output
    ans
}

   




