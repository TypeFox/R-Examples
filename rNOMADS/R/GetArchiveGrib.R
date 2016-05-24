ArchiveGribGrab <- function(abbrev, model.date, model.run, pred, local.dir = ".", file.name = "fcst.grb",
    tidy = FALSE, verbose = TRUE, download.method = NULL, file.type = "grib2") {
    #Get archived grib file
    #INPUTS
    #    ABBREV - Model abbreviation returned by NOMADSArchiveList
    #    MODEL.DATE - The year. month, and day of the model run, in YYYYMMDD format
    #    MODEL.RUN - Which hour the model was run (i.e. 00, 06, 12, 18 for GFS)
    #    PRED - Which prediction to get (analysis is 00)
    #    LOCAL.DIR is the directory to save the files in
    #    FILE.NAME is the directory path and file name to save the grib file on disk, defaults to "fcst.grb" in current directory
    #    TIDY asks whether to delete all grib files in the directory specified in FILE.NAME, default FALSE.
    #    This is useful to clear out previous model runs.
    #    It looks for all files named '.grb' and removes them.
    #    VERBOSE gives a blow by blow account of the download. Default TRUE.
    #    DOWNLOAD.METHOD allows the user to set the download method used by download.file
    #    FILE.TYPE specifies which file type to try and get - there may be only one available.
    #       grib1 for GRIB1, grib2 for GRIB2
    #OUTPUTS
    #    GRIB.INFO contains information about the downloaded file
    #        $LOCAL.DIR is the directory where the grib file is saved
    #        $FILE.NAME is the local file name where the data is stored
    #        $URL is the URL from which the file was downloaded

   if(tidy) {
        unlink(list.files(local.dir, pattern = "*\\.grb$"))
   }

    model.date <- as.numeric(strsplit(as.character(model.date), split = "")[[1]])

    if(length(model.date) != 8) {
        stop("MODEL.DATE should be in YYYYMMDD format!")
    }

    #Get specified file format
    if(file.type == "grib1") {
        suffix <- ".grb"
    } else if (file.type == "grib2") {
       suffix <- ".grb2"
    } else {
        stop(paste0("Did not recognize file type \"", file.type, "\""))
    }

    #Get model info and set up URL to archive
    model.url <- NOMADSArchiveList("grib", abbrev=abbrev)$url
    download.url <- paste0(model.url, paste(model.date[1:6], collapse = ""), "/", paste(model.date, collapse = ""), "/") 
    file.part <- paste0(paste(model.date, collapse = ""), "_", 
    sprintf("%02.f", as.numeric(model.run)), 
    "00_", sprintf("%03.f", as.numeric(pred)), suffix)

    #Find out which grib files are in the archive
    link.list <- unique(LinkExtractor(download.url))

    #Check if the requested file is where we think it is
    if(sum(grepl(paste0(".*", file.part, "$"), link.list)) < 1) {
        stop(paste("The requested data file ending in", file.part, "does not appear to be in the archive. 
            Try opening", download.url, "in a web browser to verify that it's missing."))
    }

    #Set up URL to file
    grb.urls <- paste0(download.url, link.list[grepl(paste0(".*", file.part, "$"), link.list)])

    #Download the file
    if(length(grb.urls) > 1) {
        warning("Two files were found:  ", paste(link.list[grepl(paste0(".*", file.part, "$"), link.list)], collapse = " "), ".  Both will be downloaded.")
        c <- 1 
    } else {
        c <- "" 
    }

    use.curl <- FALSE #May use this as an option in the future
    file.names <- NULL
    for(grb.url in grb.urls) {
        file.name.tmp <-  paste0(file.name, c)
       if(is.null(download.method) & !use.curl) {#Let R decide how to download the file 
           download.file(grb.url, paste(local.dir,file.name.tmp, sep = "/"), mode = "wb", quiet = !verbose)
        }
        if(!is.null(download.method) & !use.curl) { #Download using specific method
            download.file(grb.url, paste(local.dir,file.name.tmp, sep = "/"), download.method, mode = "wb", quiet = !verbose)
        }
        file.names <- append(file.names, file.name.tmp)
        if(c != "") {
            c <- c + 1
        }
    }
    grib.info <- list(local.dir = normalizePath(local.dir), file.name = file.names, url = grb.urls)
    return(grib.info)
} 

CheckNOMADSArchive <- function(abbrev, model.date = NULL) {
    #Determine what data are available in the archive for a given model.
    #If MODEL.DATE is NULL, find out what dates have directories (and probably data).
    #IF MODEL.DATE is a valid date, find out what files are available for that specific directory.
    #INPUTS
    #    ABBREV is the model abbreviation that rNOMADS uses to figure out which model you want.
    #    MODEL.DATE is a specific date in YYYYMMMDD format - if this is not NULL, figure out what model files are available for that date
    #OUTPUTS
    #    AVAILABLE.MODELS
    #        $DATE - What date the file is for, in YYYYMMDD format
    #        $MODEL.RUN - At which hour the model was run (GMT)
    #        $PRED - What predictions are available
    #        $FILE.NAME - The grib file with model data for the date, model run, and prediction

    model.url <- NOMADSArchiveList("grib", abbrev=abbrev)$url
    model.list <- c()
    if(is.null(model.date)) { #Check all available dates
        #Find out which months are available
        month.list <- grep("\\d{6}/", LinkExtractor(model.url), value = TRUE)
        for(month in month.list) {
            date.list <- grep("\\d{8}/", LinkExtractor(paste0(model.url, month)), value = TRUE)
            for(day in date.list) {
                model.list <- append(model.list, grep("grb\\d?$", LinkExtractor(paste0(model.url, month, day)), value = TRUE))
             }
         }
     } else {
            model.date <- as.numeric(strsplit(as.character(model.date), split = "")[[1]])
            model.list <- append(
                model.list, 
                grep("grb\\d?$", 
                LinkExtractor(paste0(model.url, paste(model.date[1:6], collapse = ""), "/", paste(model.date, collapse = ""), "/")),
                value = TRUE))       
     }

     model.list <- as.vector(model.list)

     available.models <- list(
         date = stringr::str_extract(model.list, "\\d{8}"),
         model.run = stringr::str_replace_all(stringr::str_extract(model.list, "_\\d{4}_"), "_", ""),
         pred = stringr::str_replace(stringr::str_replace(stringr::str_extract(model.list, "_\\d{3}\\."), "\\.", ""), "_", ""),
         file.name = model.list)
     invisible(available.models)
        
}
