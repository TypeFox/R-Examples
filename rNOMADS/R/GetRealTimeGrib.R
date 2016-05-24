CrawlModels <- function(abbrev = NULL, url = NULL, depth = NULL, verbose = TRUE) {
   #A simple web crawler that looks at the specified model directory online and gets information on all runs of the specified model.
   #See the NOMADSRealTimeList function for available models.
   #Alternatively, pass CrawlModels a URL to get a model that I have not included yet.
   #INPUTS
   #    ABBREV - Model abbreviation as defined in NOMADSRealTimeList().  
   #        If NULL, use the url you provided, if you did not provide one, throw error.
   #    URL - Use your own URL and attempt to get model data from it.  
   #        This is in case NOMADS updates its system before I have a chance to update rNOMADS
   #    DEPTH - How many links to return; set this to 1 if you only want the latest model (this will speed things up significantly)
   #    VERBOSE - Print each link you find as you find it
   #OUTPUTS
   #    URLS.OUT is a list of available models from the given ABBREV or URL

   if(is.null(url) & is.null(abbrev)) {
       stop("No models specified.")
   }
   
   if(is.null(url)) {
       model.info <- NOMADSRealTimeList("grib", abbrev=abbrev) 
       url <- model.info$url[1]
   }   

   urls.out <- unlist(WebCrawler(url, depth = depth, verbose = verbose), recursive = TRUE, use.names = FALSE) 
}

GribGrab <- function(model.url, pred, levels, variables, local.dir = ".", file.name = "fcst.grb", 
    model.domain = NULL, tidy = FALSE, verbose = TRUE, check.url = TRUE, download.method = NULL)
{
    #Get grib file from the GFS forecast repository
    #INPUTS
    #    MODEL.URL is the URL of the model, as one of the elements returned by CrawlModels
    #    PRED is the exact model run you want, generally from ParseModelPage
    #    LEVELS is the vertical region to return data for,  as vector, generally from ParseModelPage
    #    VARIABLES is the data to return, as vector, generally from ParseModelPage
    #    LOCAL.DIR is the directory to save the files in
    #    FILE.NAME is the directory path and file name to save the grib file on disk, defaults to "fcst.grb" in current directory
    #    MODEL.DOMAIN is a vector of latitudes and longitudes that specify the area to return a forecast for
    #    This is a rectangle with elements: west longitude, east longitude, north latitude, south latitude
    #    Defaults to entire planet
    #    TIDY asks whether to delete all grib files in the directory specified in FILE.NAME, default FALSE.
    #    This is useful to clear out previous model runs.
    #    It looks for all files named '.grb' and removes them.
    #    VERBOSE gives a blow by blow account of the download. Default TRUE.
    #    CHECK.URL verifies that MODEL.URL is real and contains variable and level data using ParseModelPage.
    #        If it finds that one or more levels and/or variables are missing from the URL, it throws a warning.
    #        If there is no level or variable information, or PRED is not in the list, it throws an error and exits
    #    DOWNLOAD.METHOD allows the user to set how download.file accesses the Grib file (e. g. using "internal", "wget", "curl", or "lynx"), defaults to NULL (let R decide)
    #OUTPUTS
    #    GRIB.INFO contains information about the downloaded file
    #        $LOCAL.DIR is the directory where the grib file is saved
    #        $FILE.NAME is the local file name where the data is stored
    #        $URL is the url used to retrieve the data 

   if(check.url) {
       test.scan <- ParseModelPage(model.url)
       if(length(test.scan$pred) == 0 & length(test.scan$levels) == 0 & length(test.scan$levels) == 0) {
           stop("There is no prediction, level, or variable data present in the specified model URL.
               Perhaps the model has not been loaded yet; try using an earlier model run.")
       }
       if(!(pred %in% test.scan$pred)) {
           stop(paste0("The requested prediction \"", pred, "\" was not found on the model download website."))
       }
       for(lvl in levels) {
           if(!(stringr::str_replace_all(lvl, " ", "_")  %in% test.scan$levels)) {
               warning(paste0("Requested level \"", lvl, "\" was not found on the model download website."))
           }
       }
       for(var in variables) {
           if(!(var %in% test.scan$var)) {
               warning(paste0("Requested variable \"", var, "\" was not found on the model download website."))
           }
       }
   }
   if(tidy) {
        unlink(list.files(local.dir, pattern = "*\\.grb$"))
   }


   model.str <- strsplit(model.url, "?dir=")[[1]]
   if(length(levels) > 0 & !is.null(levels)) {
        levels.str <- paste0("&lev_", paste(gsub(" ", "_", SanitizeURL(levels)), collapse = "=on&lev_"), "=on")
   } else {
       levels.str <- ""
   }
   variables.str <- paste(SanitizeURL(variables), collapse = "=on&var_")

   if(!is.null(model.domain)) {
       subregion.str <- paste("=on&subregion=",
       "&leftlon=", model.domain[1],
       "&rightlon=", model.domain[2],
       "&toplat=", model.domain[3],
       "&bottomlat=", model.domain[4],
       "&", sep = "")
    } else {
       subregion.str <- "=on&" 
    }

   grb.url <- paste0(paste0(model.str[1], "file=", pred),
       levels.str,
       "&var_",
       variables.str,
       subregion.str,
       paste0("dir=", model.str[2]))

   #now write download logic
   use.curl <- FALSE #May use this as an option in the future 
   if(is.null(download.method) & !use.curl) {#Let R decide how to download the file 
       download.file(grb.url, paste(local.dir,file.name, sep = "/"), mode = "wb", quiet = !verbose)
   }
   if(!is.null(download.method) & !use.curl) { #Download using specific method
       download.file(grb.url, paste(local.dir,file.name, sep = "/"), download.method, mode = "wb", quiet = !verbose) 
   }

   grib.info <- list(local.dir = normalizePath(local.dir), file.name = file.name, url = grb.url) 
   return(grib.info)
}

ParseModelPage <- function(model.url) {
#    This function determines the available predictions, levels, and variables for a given model page.
#    It returns a list of these predictions, levels, and variables so that a call to the model can be constructed.
#    INPUTS
#        MODEL.URL is one of the model pages returned by CrawlModels
#    OUTPUTS
#        MODEL.PARAMETERS - a list with elements
#            MODEL.PARAMETERS$PRED - Individual "predictions" - i. e. individual outputs for each model instance
#            MODEL.PARAMETERS$LEVELS - the model levels
#            MODEL.PARAMETERS$VARIABLES - the types of data provided by the models

    html.code <- scrapeR::scrape(model.url, parse = FALSE)
    model.parameters <- list()
    model.parameters$pred <- gsub("option value=\"", "", stringr::str_extract_all(html.code[[1]], "option value=\"[^\"]+")[[1]])
    model.parameters$levels <- gsub("lev_", "", stringr::str_extract_all(html.code[[1]], "lev_[^\"]+")[[1]], fixed = TRUE)
    model.parameters$variables <- gsub("var_", "", stringr::str_extract_all(html.code[[1]], "var_[^\"]+")[[1]], fixed = TRUE)
     
    return(model.parameters)
}

SanitizeURL <- function(bad.strs) {
#    This function replaces illegal characters in levels and variables prior to constructing the model retrieval URL.
#    INPUTS
#        BAD.STRS - A vector of strings, possibly with illegal URL characters
#    OUTPUTS
#        GOOD.STRS - A vector of strings with illegal characters replaced by URL codes

   good.strs <- stringr::str_replace_all(bad.strs,  "\\^", "%5E")
   good.strs <- stringr::str_replace_all(good.strs,  "\\\\\\(", "%5C%28")
   good.strs <- stringr::str_replace_all(good.strs,  "\\\\\\)", "%5C%29")
   good.strs <- stringr::str_replace_all(good.strs,  "\\\\", "%5C")
   good.strs <- stringr::str_replace_all(good.strs, "=","%3D")
   good.strs <- stringr::str_replace_all(good.strs, "/", "%2F")

   return(good.strs)
}

WebCrawler <- function(url, depth = NULL, verbose = TRUE) {
#    This function recursively searches for links in the given url and follows every single link.
#    It returns a list of the final (dead end) URLs.
#    Many thanks to users David F and Adam Smith on stackoverflow for the link parser:
#    http://stackoverflow.com/questions/3746256/extract-links-from-webpage-using-r/3746290#3746290
#    INPUTS
#        URL is the url to start looking in
#    OUTPUTS
#        URLS.OUT are the URLs at the end of the road

    links <- LinkExtractor(url)
    if(length(links) == 0) {
        if(verbose) {
            print(url)
        }
        return(url)
    } else {
        urls.out <- vector("list", length = length(links))
        for(link in links) {
           if(!is.null(depth)) {
               if(length(unlist(urls.out)) >= depth) {
                   break
               }
            }
           urls.out[[link]] <- WebCrawler(link, depth = depth, verbose = verbose)
        }
        return(urls.out)
    }
}
