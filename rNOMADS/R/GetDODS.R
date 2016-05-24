#use the GrADS-DODS capability of NOMADS to get ascii data

GetDODSDates <- function(abbrev, archive = FALSE, request.sleep = 1) {
    #Checks the GrADS data server to see what dates and model subsets are available for model specified by ABBREV.
    #INPUTS
    #    ABBREV - Model abbreviation
    #    ARCHIVE - If you're looking in the model archives (TRUE) or the real time NOMADS system (FALSE)
    #    REQUEST.SLEEP - Sometimes hammering the NOMADS server with a zillion HTTP requests is not a good idea.
    #    REQUEST.SLEEP pauses X seconds between requests to prevent timeouts.
    #OUTPUTS
    #    AVAILABLE.DATES - A list of model URLS and dates
    #        $ABBREV - Model abbreviation
    #        $DATE - Model run date, YYYYMMDD
    #        $URL - Model url

    date.pattern <- "[1-2]\\d{3}[0-1]\\d{1}[0-3]\\d{1}$"
    
    if(!archive) {
        top.url <- unique(NOMADSRealTimeList("dods", abbrev)$url)
    } else {
        if(grepl("anl$", abbrev)) {
            stop(paste("Archived analysis models are not stored by date.",
                  "To find out the available model runs, get the model URL from NOMADSArchiveList and pass it to GetDODSModelRuns.",
                  "Then use GETDODSModelRunInfo to find out what coverage the models have, and use DODSGrab to grab the data, where the first argument is the URL from NOMADSArchiveList",
                  "and the second is the model run from GetDODSModelRuns.")
                  )
        } else {
            top.url <- unique(NOMADSArchiveList("dods",abbrev)$url)
            if(top.url == "NONE") {
                stop("The archived model you requested is not available on the NCEP DODS system.")
            }
        }
    } 
     
   if(!RCurl::url.exists(top.url)) {
       stop(paste0("The specified URL does not exist!  Make sure your model information is correct.  It is also possible the NOMADS server is down.\n",
          "Details:  Attempted to access ", top.url, " but did not succeed..."))
   }

    top.links <-  LinkExtractor(top.url)

    #List entries in html.tmp, see if they are dates or not
    date.links.lind <- grepl(date.pattern, as.character(top.links))
    if(sum(date.links.lind) == 0) { #If there do not appear to be dates here, go down one more directory level
        urls.tmp <- c()
        for(k in seq(4, length(top.links) - 3)) {
               if(!RCurl::url.exists(top.links[k])) {
                   stop(paste0("The specified URL does not exist!  Make sure your model information is correct.  It is also possible the NOMADS server is down.\n",
                       "Details:  Attempted to access ", top.links[k], " but did not succeed..."))
                }
            links.low <- LinkExtractor(top.links[k])
            urls.tmp <- append(urls.tmp, links.low[grepl(date.pattern, links.low)])
            if(request.sleep > 0) {
                Sys.sleep(request.sleep)
            }
        }
    } else {
        urls.tmp <- top.links[date.links.lind]
    }

    return(list(model = abbrev, date = stringr::str_extract(urls.tmp, date.pattern), url = urls.tmp)) 
}

GetDODSModelRuns <- function(model.url) {
   #Given a URL of a certain model date, determine which model runs are available for that date.
   #This is useful to check to see if a certain model run is on the website before attempting to get information or data from it.
   #The model url probably comes from GetDODSDates.
   #INPUTS 
   #    MODEL.URL - A URL pointing to the DODS model page for a certain date; get this URL from GetDODSDates.
   #OUTPUTS
   #    AVAILABLE.MODEL.RUNS - Which model runs are available for that date
   #        $MODEL.RUN - The model run
   #        $MODEL.RUN.INFO - Info about the model run, hence the name

      if(!RCurl::url.exists(model.url)) {
       stop(paste0("The specified URL does not exist!  Make sure your model information is correct.  It is also possible the NOMADS server is down.\n",
          "Details:  Attempted to access ", model.url, " but did not succeed..."))
   }

   html.tmp <- XML::htmlParse(model.url)
   model.runs <- XML::xpathSApply(html.tmp, '//b', XML::xmlValue) 
   XML::free(html.tmp)
   html.txt <- readLines(model.url)

   model.run.str <- NULL
   model.info.str <- NULL
   for(k in seq_len(length(model.runs))) {
       model.run.tmp <-  paste0(stringr::str_extract_all(stringr::str_replace(model.runs[k], "^\\d+:", ""), "[a-zA-Z0-9_.]")[[1]], collapse = "")
       model.info.tmp <- stringr::str_replace(html.txt[grepl(paste0(model.run.tmp, ":"), html.txt)], "</b>&nbsp;", " ")
       model.run.str <- append(model.run.str, model.run.tmp)
       model.info.str <- append(model.info.str, model.info.tmp)
   } 

   return(list(model.run = model.run.str, model.run.info = model.info.str))
}

GetDODSModelRunInfo <- function(model.url, model.run) {
   #Get description of the model run,  documentation (if present), longitude and latitude covered, time span, levels (if present), and variable list
   #INPUTS
   #    MODEL.URL is a URL pointing to a certain model date, probably from GetDODSDates.
   #    MODEL.RUN is a specified model run, probably from GetDODSModelRuns.
   #OUTPUTS
   #    MODEL.PARAMETERS - List of variables and model coverage information

   info.url <- paste0(model.url, "/", model.run, ".info")

   if(!RCurl::url.exists(info.url)) {
       stop(paste0("The specified URL does not exist!  Make sure your model information is correct.  It is also possible the NOMADS server is down.\n",
          "Details:  Attempted to access ", info.url, " but did not succeed..."))
   }
            
   info.table <- XML::readHTMLTable(info.url)[[2]]
   info.arr <- cbind(
       as.vector(info.table[,1]),
       as.vector(info.table[,2]),
       as.vector(info.table[,3]))

   info.arr[which(is.na(info.arr), arr.ind=TRUE)] <- ""
   #Get rid of bad characters
   model.info <- stringr::str_replace_all(apply(info.arr, 1, paste, collapse = " "), "\\N{LATIN CAPITAL LETTER A WITH CIRCUMFLEX}", "")
   return(model.info)
}

DODSGrab <- function(model.url, model.run, variables, time, lon, lat, levels = NULL, display.url = TRUE, verbose = FALSE, request.sleep = 1) {
   #Get data from DODS.  Note that this is slower than GribGrab but will work on all operating systems.
   #The output of this function will be the same as the output of ReadGrib in order to maintain consistency across rNOMADS.
   #ALL INDICES START FROM ZERO 
   #INPUTS
   #    MODEL.URL is a URL pointing to a certain model date, probably from GetDODSDates. 
   #    MODEL.RUN is a specified model run, probably from GetDODSModelRuns.
   #    VARIABLES is a list of variables from the list returned by GetDODSModelRunInfo
   #    TIME is an **index list** of times per info from GetDODSModelRunInfo, "c(x,y)"
   #    LON is an **index list** of longitudes per info from GetDODSModelRunInfo "c(x,y)"
   #    LAT is an **index list** of latitudes per info from GetDODSModelRunInfo "c(x,y)"
   #    LEVELS is an **index list** of levels per info from GetDODSModelRunInfo "c(x,y)"
   #         if not NULL, try to request the variable at a certain level.  Will fail if the variable does not have associated levels.
   #    DISPLAY.URL asks whether to display the URL request for debugging purposes.
   #        You can paste it into your browser to check to make sure things are working correctly
   #    VERBOSE gives a very talkative description of the download process
   #   REQUEST.SLEEP says how many seconds to pause between data requests to avoid server timeouts
   #OUTPUTS
   #    MODEL.DATA - the model as an array, with columns for the model run date (when the model was run)
   #       the forecast (when the model was for), the variable (what kind of data), the level (where in the atmosphere or the Earth, vertically)
   #       the longitude, the latitude, and the value of the variable.

   prev.digits <- options("digits")
   options("digits" = 8)

   model.data <- list(
       model.run.date = NULL,
       forecast.date  = NULL,
       variables      = NULL,
       levels         = NULL,
       lon            = NULL,
       lat            = NULL,
       value          = NULL,
       request.url    = NULL)
       
   for(variable in variables) { 
       preamble <- paste0(model.url, "/", model.run, ".ascii?", variable)
       time.str <- paste0("[", paste0(time, collapse = ":"), "]")
    
       l.ind <- !is.null(levels)
    
       if(l.ind) {
           level.str <- paste0("[", paste0(levels, collapse = ":"), "]")
       } else {
           level.str <- ""
       }
       lat.str <- paste0("[", paste0(lat, collapse = ":"), "]")
       lon.str <- paste0("[", paste0(lon, collapse = ":"), "]")
      
       data.url <- paste0(preamble, time.str, level.str, lat.str, lon.str)  
       
       if(display.url) {
           print(data.url)
       }
    
       #RCurl needs to be loaded for this to work I think
       #data.txt <- readLines(data.url)
       data.txt.raw <- RCurl::getURL(data.url, .opts = list(verbose = verbose)) #Read in data
       if(grepl("[eE][rR][rR][oO][rR]", data.txt.raw)) {
           warning(paste0("There may have been an error retrieving data from the NOMADS server.  HTML text is as follows\n", data.txt.raw
           ))
       }
       data.txt <- unlist(strsplit(data.txt.raw, split = "\\n"))
       lats <- as.numeric(unlist(strsplit(data.txt[grep("^lat,", data.txt) + 1], split = ",")))
       lons <- as.numeric(unlist(strsplit(data.txt[grep("^lon,", data.txt) + 1], split = ",")))
     
       if(l.ind) {
           levels.out <- as.numeric(unlist(strsplit(data.txt[grep("^lev,", data.txt) + 1], split = ",")))
       }
       t.ind <- grep("^time,", data.txt)

       prev.digits <- options("digits")
       options("digits" = 15)
       num.times <- as.numeric(unlist(strsplit(data.txt[t.ind + 1], split = ",")))
       options("digits" = prev.digits$digits)
       times <- as.POSIXlt(as.Date(num.times - 2, origin = "1-1-1"), tz = "GMT") + 3600 * 24 * (num.times - floor(num.times))

       #Extract data values 
       
       val.dim <- unlist(stringr::str_extract_all(unlist(strsplit(data.txt[1], ","))[2], "\\d"))
       val.txt <- data.txt[2:(t.ind - 1)]    
       val.txt <- val.txt[val.txt !=""] 
       val.txt <- stringr::str_replace_all(val.txt, "\\]\\[", ",")
       val.txt <- stringr::str_replace_all(val.txt, c("\\]|\\["), "")
    
       model.run.date <- paste0(stringr::str_extract(model.url, "[1-2]\\d{3}[0-1]\\d{1}[0-3]\\d{1}$"), model.run)
       
       row.num <- (as.numeric(val.dim[2]) + l.ind) * length(val.txt)
       model.data.tmp <- array("", dim = c(row.num, 7))
    
       r.start <- 3 + l.ind #What row to start at
       for(k in seq_len(length(val.txt))) {
           val.tmp <- sapply(strsplit(val.txt[k], split = ","), as.numeric)
           r.end <- length(val.tmp)
           get.rows <- r.end - r.start + 1
           model.data$model.run.date <- append(model.data$model.run.date, rep(model.run.date, get.rows))
           model.data$forecast.date  <- append(model.data$forecast.date, rep(times[val.tmp[1] + 1], get.rows))
           model.data$variables      <- append(model.data$variables, rep(variable, get.rows))
           if(l.ind) {
               model.data$levels <- append(model.data$levels, rep(levels.out[val.tmp[2] + 1], get.rows))
           }
           model.data$lon            <- append(model.data$lon, lons) 
           model.data$lat            <- append(model.data$lat, rep(lats[val.tmp[2 + l.ind] + 1], get.rows))
           model.data$value          <- append(model.data$value, val.tmp[r.start:r.end])
       }
       
       if(!l.ind) {
           model.data$levels <- rep("level not defined", length(model.data$value))
       }
    
       model.data$request.url <- append(model.data$request.url, data.url)
    
       if(length(variables) > 1) {
           Sys.sleep(request.sleep)
       }
}
   return(model.data)
}

