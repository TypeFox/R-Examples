##########################
##########################
#various pems imports
##########################
##########################


#additional 

#work in progress


#############################
#############################
##importRoyalTek2PEMS
#############################
#############################

#quick importer for Leed City Council GPS
#started 05-11-2011
#kr v0.0.2 (09-11-2011)

#three functions
###################
#importRoyalTek2PEMS
#importRoyalTekTXT
#importRoyalTekNMEA
#

#currently only exporting the wrapper

###################
#to do
###################
#unit assignment
##
#lat, long correction
#for north east, west, south
#

###################
#notes
###################
#

###################
#importRoyalTek2PEMS
###################

#parent function
#(TXT/NMEA switch)

importRoyalTek2PEMS <- function(file.name = file.choose(), file.type = c("special", "txt", "nmea"),
                                vbox = "RoyalTEk", history = NULL, constants = NULL, ... ){
   
    #setup
    this.call <- match.call()
    fun.name <- "importRoyalTek2PEMS"

    file.type <- checkOption(file.type[1], formals(importRoyalTek2PEMS)$file.type, 
                              "file.type", "known import methods", 
                              fun.name = fun.name)

    #set file.type if null, guessimate format
    if(file.type == "special"){
        temp <- tolower(substr(file.name, nchar(file.name)-3, nchar(file.name)))
        if(temp == ".txt") 
            file.type <- "txt"
        if(temp == "nmea")
            file.type <- "nmea"
    }

    #import is type valid
    ans <- NULL
    units <- NA
    if(is.character(file.type)){
        if(file.type == "txt")
            ans <- importRoyalTekTXT(file.name = file.name, ...)
        if(file.type == "nmea")
            ans <- importRoyalTekNMEA(file.name = file.name, ...)
    }

    #stop is not valid
    if(is.null(ans))
        stop(paste("\t In ", fun.name,"(...) selected file type not recognised", sep=""),
             paste("\n\t [suggest setting file.type in call if known]", sep=""), 
             call. = FALSE, domain = NA)

    output <- makePEMS(x = ans, units = NULL, constants = constants, history = history, ...)
    output$vbox <- vbox
    
    #reset history?
    output$history[[length(output$history)]] <- this.call 

    #return output
    invisible(output)
}


####################
#importRoyalTekTXT
####################

importRoyalTekTXT <- function(
         file.name = file.choose(), n = NULL, to.lower = TRUE,
         fields = c("Record", "Event Type", "Year", "Month", "Day", 
                    "Hour", "Minute", "Second", "Latitude", "Longitude",    
                    "Altitude", "PDOP", "HDOP", "Satellite No",        
                    "Speed", "Direction" ),
         info = c("Datalogs", "Date", "Time", "Device ID", 
                  "About total Numbers"),
         exclude = c(": ", "\\(KMs/hr\\)"),         
         ...
){

    #set up
    if(is.null(n))
        n <- -1L

    #import
    source <- readLines(file.name, n = n, ...)

    #recover columns and strip field names
    a <- sapply(fields, function(x){ 
                            temp <- source[grep(x, source)]
                            gsub(x, "", temp)})

    #strip excludes
    for(i in exclude)
        a <- gsub(i, "", a)

    #make data frame
    a <- data.frame(a, stringsAsFactors = FALSE)

    #get date
    date <- paste(a$Year, a$Month, a$Day, sep="-")
    date <- paste(date, paste(a$Hour, a$Minute, a$Second, sep="-"), sep=" ")

    #make basic data numeric
    a <- as.data.frame(cbind(date = date, a), stringsAsFactors = FALSE)

    #set date format
    a$date <- as.POSIXct(strptime(a$date, format = "%Y-%m-%d %H-%M-%OS", "GMT"))

    a[,2:ncol(a)] <- apply(a[,2:ncol(a)], 2, as.numeric)

    #name to lower?
    if(to.lower)
        names(a)<-tolower(names(a))

    #get info
    b <- sapply(info, function(x) source[grep(x, source)])
    comment(a) <- paste(b)

    #output
    invisible(a)
}


#####################
#importRoyalTekNMEA
#####################

#quick importer for Leed City Council GPS
#NMEA handler

importRoyalTekNMEA <- function(
         file.name = file.choose(), n = NULL, to.lower = TRUE,
         fields = NULL, filter.method = NULL,
         filter = "\\$GPGGA",         
         ...
){

    #set up
    if(is.null(n))
        n <- -1L

    #import
    source <- readLines(file.name, n = n, ...)

    #filter

    if(filter == "\\$GPGGA"){
        if(is.null(filter.method)) 
            filter.method <- c("del", "time", "lat", "f", "long", "f", "n", "n", "n", "n", "c", "c", "c", "c", "c")
        if(is.null(fields))
            fields <- c("time", "latitude", "North", "longitude", "West", "X", "X", "X", "X", "X", "X", "X", "X", "X")
    }

    #stop if invalid filter   
    if(is.null(filter.method) | is.null(fields))
        stop("INVALID FILTER/FIELDS")

    ans <- source[grep(filter[1], source)]

    #make dataframe for filter
#comma delim currently hard coded
    ans <- data.frame(do.call(rbind, strsplit(ans, ",")), stringsAsFactors = FALSE)

    #drop deleted cases 
    ans <- ans[filter.method != "del"]
    filter.method <- filter.method[filter.method != "del"]

    #rename
    names(ans) <- make.names(fields, unique = TRUE)
    if(to.lower)
        names(ans) <- tolower(names(ans))

#finish below
#validate below
#rationalise below

    #get list of times
#work to do here
#need date for full time stamp
    temp <- names(ans)[filter.method=="time"]
    for(i in temp){
        ans[, i] <- as.numeric(ans[, i])
    }

    #get list of lat, long
#needs validating
    temp <- names(ans)[filter.method=="lat" | filter.method=="long"]
    for(i in temp){
        unit <- as.numeric(ans[, i])
        ans[, i] <- floor(unit / 100)
        unit <- unit - (ans[, i] * 100)
#        ans[, i] <- ans[, i] + (floor(unit) / 60)
#        ans[, i] <- ans[, i] + (unit - floor(unit)) / 60
        ans[, i] <- ans[, i] + ((unit) / 60)

    }

    #get list of numerics
    temp <- names(ans)[filter.method=="n"]
    for(i in temp){
        ans[, i] <- as.numeric(ans[, i])
    }

    #get list of characters
    temp <- names(ans)[filter.method=="c"]
    for(i in temp){
        ans[, i] <- as.character(ans[, i])
    }

    #get list of factors
    temp <- names(ans)[filter.method=="f"]
    for(i in temp){
        ans[, i] <- as.factor(ans[, i])
    }

    #output
    invisible(ans)

}

