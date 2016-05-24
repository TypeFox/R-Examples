FDSNGetEvents <- function(base.url, parameters, verbose = TRUE) {
    #Get event infromation from teh International Federation of Digital Seismograph Networks (FDSN).
    #INPUTS
    #    BASE.URL - Where the FDSN server is located, for example http://service.iris.edu/fdsnws/event/1/
    #    PARAMETERS - Selection parameters 
    #         $NAME - Name of parameter (see available ones via FDSNParameterList("event"))
    #         $VALUE - Value of parameter
    #    VERBOSE - If TRUE, print constructed URL
    #OUTPUTS
    #    RESULT - Earthquake information returned by the query.

   if("format" %in% parameters$name) {
        if(parameters$value[parameters$name == "format"] == "text") {
            warning("The rFDSN package requires an XML download format...and it will probably crash now.  Please set the format back to xml.")
        }
    }

    if(length(parameters$value) != length(parameters$name)) {
        stop("The number of parameters and the number of parameter values are not equal!  Please check inputs.")
    }

    url.to.get <- paste0(base.url, "query?", paste(parameters$name, parameters$value, sep = "=", collapse = "&"))
    if(verbose) {
        print(url.to.get)
    }
    doc <- xmlInternalTreeParse(url.to.get)
    top <- xmlChildren(xmlRoot(doc))[1]

    event <- list(
        publicID =              NULL,
        type     =              NULL,
        descriptionType =       NULL,
        descriptionText =       NULL,
        preferredMagnitudeID  = NULL,
        preferredOriginID =     NULL,
        origin =                NULL,
        magnitude =             NULL
     ) 

   top <- xmlChildren(xmlRoot(doc))$eventParameters
   n.i <- 1
   for(k in which(names(xmlApply(top, xmlChildren)) == "event")) {
       event$publicID <- append(event$publicID, as.character(xpathApply(top[[k]], "@publicID")))
       event$type <- append(event$type, xmlValue(top[[k]]["type"]$type))
       event$descriptionType <- append(event$descriptionType, xmlValue(xmlSApply(top[[k]], xmlChildren)$description$type))
       event$descriptionText <- append(event$descriptionText, xmlValue(xmlSApply(top[[k]], xmlChildren)$description$text))
       event$preferredMagnitudeID <- append(event$preferredMagnitudeID, xmlValue(top[[k]]["preferredMagnitudeID"]$preferredMagnitudeID))
       event$preferredOriginID <- append(event$preferredOriginID, xmlValue(top[[k]]["preferredOriginID"]$preferredOriginID))
       
       origin <- list(
           publicID =                  NULL,
           iris.contributor =          NULL,
           iris.catalog =              NULL,
           iris.contributorOriginId =  NULL,
           iris.contributorEventId =   NULL,
           time =                      NULL,
           author =                    NULL,
           latitude =                  NULL,
           longitude =                 NULL,
           depth =                     NULL
       )           

       for(j in which(names(xmlApply(top[[k]], xmlChildren)) == "origin")) {
           ortmp <- xmlChildren(top[[k]])[j]$origin
           origin$publicID <- append(origin$publicID, as.character(xpathApply(ortmp, "@publicID"))) 
           origin$iris.contributor <- append(origin$iris.contributor, as.character(xpathApply(ortmp, "@iris:contributor")))
           origin$iris.catalog <- append(origin$iris.catalog, as.character(xpathApply(ortmp, "@iris:catalog")))
           origin$iris.contributorOriginId <- append(origin$iris.contributorOriginId, as.character(xpathApply(ortmp, "@iris:contributorOriginId")))
           origin$iris.contributorEventId <- append(origin$iris.contributorEventId, as.character(xpathApply(ortmp, "@iris:contributorEventId")))
           origin[["time"]] <- append(origin[["time"]], xmlValue(ortmp[["time"]][["value"]]))
           origin[["author"]] <- append(origin[["author"]], xmlValue(ortmp[["creationInfo"]][["author"]]))
           origin[["latitude"]] <- append(origin[["latitude"]], xmlValue(ortmp[["latitude"]][["value"]]))
           origin[["longitude"]] <- append(origin[["longitude"]], xmlValue(ortmp[["longitude"]][["value"]]))
           origin[["depth"]] <- append(origin[["depth"]], xmlValue(ortmp[["depth"]][["value"]]))
       }

       event$origin[[n.i]] <- origin

        magnitude <- list(
             publicID =  NULL,
             originID =  NULL,
             type =      NULL,
             magnitude = NULL,
             author =    NULL
      ) 

      for(j in which(names(xmlApply(top[[k]], xmlChildren)) == "magnitude")) {
           matmp <- xmlChildren(top[[k]])[j]$magnitude
           magnitude$publicID <- append(magnitude$publicID, as.character(xpathApply(matmp, "@publicID")))
           magnitude$originID <- append(magnitude$originID, xmlValue(matmp[["originID"]]))
           magnitude$type <- append(magnitude$type, xmlValue(matmp[["type"]]))
           magnitude$magnitude <- append(magnitude$magnitude, xmlValue(matmp[["mag"]][["value"]]))
           magnitude$author <- append(magnitude$author, xmlValue(matmp[["creationInfo"]][["author"]]))
      }

      event$magnitude[[n.i]] <- magnitude
      n.i <- n.i + 1
    }

    invisible(event)
}
FDSNGetNetworks <- function(base.url, parameters, verbose = TRUE) {
    #Get station information from the International Federation of Digital Seismograph Networks  (FDSN).   
    #INPUTS
    #    BASE.URL - Where the FDSN server is located, for example http://service.iris.edu/fdsnws/station/1/
    #    PARAMETERS - Selection parameters 
    #         $NAME - Name of parameter (see available ones via FDSNParameterList("station"))
    #         $VALUE - Value of parameter
    #    VERBOSE - If TRUE, print constructed URL
    #OUTPUTS
    #    RESULT - Network and station information returned by the query.

    #Check parameters for "text" if so, throw error

    if("format" %in% parameters$name) {
        if(parameters$value[parameters$name == "format"] == "text") {
            warning("The rFDSN package requires an XML download format...and it will probably crash now.  Please set the format back to xml.")
        }
    }

    if(length(parameters$value) != length(parameters$name)) {
        stop("The number of parameters and the number of parameter values are not equal!  Please check inputs.")
    }

    url.to.get <- paste0(base.url, "query?", paste(parameters$name, parameters$value, sep = "=", collapse = "&"))
    if(verbose) {
        print(url.to.get)
    }
    doc <- xmlInternalTreeParse(url.to.get)
    top <- xmlRoot(doc)
    network <- list(
        code =                       NULL, 
        startDate =                  NULL,
        endDate =                    NULL,
        restrictedStatus =           NULL,             
        Description =                NULL,
        TotalNumberStations =        NULL,
        SelectedNumberStations =     NULL,
        Stations =                   NULL)
    n.i <- 1
    for(k in which(names(xmlApply(top, xmlChildren)) == "Network")) {
        #Get network attributes
        for(name in names(network)[1:4]) {
            network[[name]] <- append(network[[name]], as.character(xpathApply(top[[k]], paste0("@", name))))
        }

        #Get top level network children
        netinfo <- xmlSApply(top[[k]], xmlValue)
        for(name in names(netinfo)) {
            if(name != "Station") {
                network[[name]] <- append(network[[name]], netinfo[[name]]) 
            }
        } 
        
        #Get station information
        station <- list(
            code =                   NULL,
            startDate =              NULL,
            endDate =                NULL,
            restrictedStatus =       NULL,
            Latitude =               NULL,
            Longitude =              NULL,
            Elevation =              NULL,
            Site =                   NULL,
            CreationDate =           NULL,
            TotalNumberChannels =    NULL,
            SelectedNumberChannels = NULL)

        for(j in which(names(xmlChildren(top[[k]])) == "Station")) {
            sta <- xmlChildren(top[[k]])[j]$Station
    
            for(name in names(station)[1:4]) {
                station[[name]] <- append(station[[name]], as.character(xpathApply(sta, paste0("@", name))))
            }
    
            stainfo <- xmlSApply(sta, xmlValue)         
            for(name in names(station)) {
                if(name %in% names(xmlChildren(sta))) { 
                   station[[name]] <- append(station[[name]], stainfo[[name]])
                }
            }
      }
      network$Stations[[n.i]] <- station
      n.i <- n.i + 1
   }

   invisible(network)
}

FDSNGetTimeSeries <- function(base.url, parameters, save.file = "result.mseed", save.dir = ".", verbose = TRUE) {
    #Download requested time series from FDSN stations in miniseed format
    #INPUTS
    #    BASE.URL - Where the FDSN server is located, for example http://service.iris.edu/fdsnws/dataselect/1/
    #    PARAMETERS - Selection parameters 
    #         $NAME - Name of parameter (see available ones via FDSNParameterList("dataselect"))
    #         $VALUE - Value of parameter
    #    SAVE.FILE - Name of the miniseed file to be saved
    #    SAVE.DIR - Where to save downloaded files, defaults to current directory
    #    VERBOSE - If TRUE, print constructed URL
    #OUTPUTS
    #    FILE.NAME - File name and path of downloaded file


    if(length(parameters$value) != length(parameters$name)) {
        stop("The number of parameters and the number of parameter values are not equal!  Please check inputs.")
    }

    if(!file.exists(save.dir)) {
        stop(paste("Directory", save.dir, "does not exist!"))
    }

    url.to.get <- paste0(base.url, "query?", paste(parameters$name, parameters$value, sep = "=", collapse = "&"))
   
    options(show.error.messages = FALSE)
    err <- try(download.file(url.to.get, paste0(save.dir, "/", save.file), quiet = !verbose))
    options(show.error.messages = TRUE)
   
    #Don't come to a screeching halt if the file is not there 
    if(class(err) != "try-error") {
         return(paste0(save.dir, "/", save.file))
    } else {
         warning(err)
         return(NULL)
    }
}


