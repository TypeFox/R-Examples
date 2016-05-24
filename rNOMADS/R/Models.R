#Descriptions of real time and archived models
NOMADSRealTimeList <- function(url.type, abbrev = NULL) {
    #Returns a list of model abbreviations for real time models, a short description, and URL for each model offered by the NOMADS server
    #If a specific model abbreviation is requested, the abbreviation is checked against the model list.
    #If a match is found, information is returned about that model; otherwise an error occurs
    #
    #
    #A big shout out to user hrbrmstr on Stack Overflow for providing the table parsing code in this function
    #http://stackoverflow.com/questions/27592575/dropped-rows-using-readhtmltable-in-r
    #http://stackoverflow.com/users/1457051/hrbrmstr
    #
    #INPUTS
    #    URL.TYPE determines which URL to return: one for downloading GRIB files (grib) or one for downloading dods data via DODS (dods)
    #    ABBREV is the model abbreviation that rNOMADS uses to figure out which model you want.
    #        if NULL, returns information on all models
    #OUTPUTS
    #    MODEL.LIST - a list of model metadata with elements
    #        $ABBREV - the abbrevation used to call the model in rNOMADS
    #        $NAME - the name of the model
    #        $URL - the location of the model on the NOMADS website

    if (!(url.type %in% c("grib", "dods"))) {
        stop("URL type must be either \"grib\" or \"dods\"!")
    }

    base.url <- "http://nomads.ncep.noaa.gov/"
    trim <- function(x) gsub("^[[:space:]]+|[[:space:]]+$", "", x)

    doc <- xml2::read_html(base.url)
    ds <- doc %>% html_nodes(xpath="//table/descendant::th[@class='nomads'][1]/../../
                                            descendant::td[contains(., 'http')]/
                                            preceding-sibling::td[3]")
    data.set <- ds %>% html_text() %>% trim()

    grib.filter <- doc %>% html_nodes(xpath="//table/descendant::th[@class='nomads'][1]/../../
                                  descendant::td[contains(., 'http')]/preceding-sibling::td[1]") %>%
   sapply(function(x) {
     ifelse(grepl("href", as.character(x)),
           x %>% html_node("a") %>% html_attr("href"),
           NA)
    })

   http.link <- doc %>% html_nodes("a[href^='/pub/data/']") %>% html_attr("href")

   gds.alt <- doc %>% html_nodes(xpath="//table/descendant::th[@class='nomads'][1]/../../
                              descendant::td[contains(., 'http')]/following-sibling::td[1]") %>%
   sapply(function(x) {
     ifelse(grepl("href", as.character(x)),
           x %>% html_node("a") %>% html_attr("href"),
           NA)
    })

 
   grib.abbrevs <- stringr::str_replace(stringr::str_replace(basename(grib.filter), "filter_", ""), ".pl", "")
   dods.abbrevs <- basename(gds.alt)
   dods.base.url <- "http://nomads.ncep.noaa.gov:9090/dods/"
   if(is.null(abbrev)) {
       if(url.type == "grib") {
          good.abbrevs <- which(!is.na(grib.abbrevs))
          model.list <- list(abbrev = grib.abbrevs[good.abbrevs], 
               name = data.set[good.abbrevs], 
               url = paste0(base.url, 
               grib.filter[good.abbrevs]))
       } else {
          good.abbrevs <- which(!is.na(dods.abbrevs))
          model.list <- list(abbrev = dods.abbrevs[good.abbrevs], 
              name = data.set[good.abbrevs], 
              url = paste0(dods.base.url, 
              basename(gds.alt[good.abbrevs]), "/"))
       }
  } else {
      if(url.type == "grib") {
           abbrev.ind <- which(abbrev == grib.abbrevs)
           if(length(abbrev.ind) > 0) {
               model.list <- list(abbrev = grib.abbrevs[abbrev.ind], 
                   name = data.set[abbrev.ind], 
                   url = paste0(base.url, grib.filter[abbrev.ind]))
           } else {
                stop(paste0("The model you searched for: \"", 
                    abbrev, 
                    "\" is not included in NOMADS real time ", 
                    url.type,
                     " model products.  Sorry!"))    
           }
      } else {
           abbrev.ind <- which(abbrev == dods.abbrevs)
           if(length(abbrev.ind) > 0) {
               model.list <- list(abbrev = dods.abbrevs[abbrev.ind], 
                   name = data.set[abbrev.ind], 
                   url = paste0(dods.base.url, basename(gds.alt[abbrev.ind]), "/"))
           } else {
                stop(paste0("The model you searched for: \"", 
                    abbrev, 
                    "\" is not included in NOMADS real time ", 
                    url.type, 
                     " model products.  Sorry!"))
           }
       }
   }
   
   return(model.list)
}
   
NOMADSArchiveList <- function(url.type, abbrev = NULL) {
    #Returns a list of model abbreviations for archived models, a short description, and URL for each model offered by the NOMADS server
    #If a specific model abbreviation is requested, the abbreviation is checked against the model list.
    #If a match is found, information is returned about that model; otherwise an error occurs
    #INPUTS
    #    URL.TYPE determines which URL to return: one for downloading GRIB files (grib) or one for downloading dods data via DODS (dods)
    #    ABBREV is the model abbreviation that rNOMADS uses to figure out which model you want.
    #        If NULL, returns information on all models
    #OUTPUTS
    #    MODEL.LIST - a list of model metadata with elements
    #        $ABBREV - the abbrevation used to call the model in rNOMADS
    #        $NAME - the name of the model
    #        $URL - the location of the model on the NOMADS website

   if (!(url.type %in% c("grib", "dods"))) {
        stop("URL type must be either \"grib\" or \"dods\"!")
    }


    abbrevs <- c(
        "ruc",
        "ruc13",
        "meso-eta-hi",
        "gfs-avn-hi",
        "gfs4",
        "rap252",
        "rap130",
        "gfsanl",
        "rucanl",
        "namanl"
       )

    names <- c(
        "Rapid Update Cycle 20 km grid",
        "Rapid Update Cycle 13 km grid",
        "North American Mesoscale, Near Real Time",
        "Global Forecast System, Near Real Time, 1 degree grid",
        "Global Forecast System, Near Real Time, 0.5 degree grid",
        "Rapid Refresh Weather Prediction System - Near Real Time, 20 km grid",
        "Rapid Refresh Weather Prediction System - Near Real Time, 13 km grid",
        "Global Forecast System, Analysis",
        "Rapid Update Cycle, Analysis",
        "North American Mesoscale, Analysis"
         )

    dods.urls <- c(
        "http://nomads.ncdc.noaa.gov/dods/NCEP_RUC/",
        "NONE",
        "http://nomads.ncdc.noaa.gov/dods/NCEP_NAM/",
        "http://nomads.ncdc.noaa.gov/dods/NCEP_GFS/",
        "NONE",
        "http://nomads.ncdc.noaa.gov/dods/NCEP_RUC/",
        "NONE",
        "http://nomads.ncdc.noaa.gov/dods/NCEP_GFS_ANALYSIS/",
        "NONE",
        "http://nomads.ncdc.noaa.gov/dods/NCEP_NAM_ANALYSIS/"        
    )

    if(!is.null(abbrev)) {
        i <- which(abbrevs == abbrev)
        if(length(i) == 0) {
            stop(paste("The model you searched for:\"", abbrev, "\"is not included in rNOMADS.  Sorry!"))
        } else {
            if(url.type == "grib") {
                 url <- paste0("http://nomads.ncdc.noaa.gov/data/", abbrev, "/")
            } else {
                url <- dods.urls[i]
            }
             return(list(abbrev = abbrev, name = names[i], url = url))
        }
    }

    if(url.type == "grib") {
        url <- paste0("http://nomads.ncdc.noaa.gov/data/", abbrevs, "/")
    } else {
        url <- dods.urls
    }

    return(list(abbrevs = abbrevs, names = names, url = url))
}
