#' Function to import daily hydrologic time series 
#' data given a USGS streamgage identification number.
#'
#' This function will import data from a WaterML1 service (current USGS 
#' hydrological data standard).  It will retrieve daily streamflow and 
#' continuous water-quality data from the USGS Daily Values Site Web 
#' Service \url{http://waterservices.usgs.gov/rest/DV-Service.html}
#' (U.S. Geological Survey, 2012d).
#'
#' @name importDVs
#' @title Imports daily USGS hydrologic times series data
#' @param staid is the USGS site identification number,  
#' usually eight digits long, but can be longer.  Users may search for 
#' surface-water sites and obtain station identification numbers using the 
#' USGS Site Web Service, 
#' \url{http://waterservices.usgs.gov/rest/Site-Service.html} (USGS, 2012e); 
#' using the National Water Information System: Mapper, 
#' \url{http://wdr.water.usgs.gov/nwisgmap/} (U.S. Geological Survey, 2012a); 
#' or using the National Water Information System: Web Interface to daily 
#' surface-water data, 
#' \url{http://waterdata.usgs.gov/nwis/dv/?referred_module=sw} 
#' (U.S. Geological Survey, 2012f).  The site identification number needs to 
#' be entered as a character, that is in quotes, because many USGS 
#' streamgage numbers begin with zero and the leading zero is necessary.
#' @param code is the USGS parameter code, a 5-digit number 
#' used in the USGS computerized data system, National Water 
#' Information System (NWIS), to uniquely identify a specific hydrologic 
#' property or constituent.  A list of paramater codes is available at 
#' \url{http://nwis.waterdata.usgs.gov/usa/nwis/pmcodes} (U.S. Geological Survey, 2012b).
#' @param stat is the USGS statistics code, a 5-digit number 
#' used in the USGS computerized data system, National Water 
#' Information System (NWIS), to uniquely identify specific statistics, such
#' as daily mean, daily maximum, and daily minimum.  The default,  
#' 00003,  is the mean daily value.  A list of statistics codes is available at 
#' \url{http://nwis.waterdata.usgs.gov/nwis/help/?read_file=stat&format=table} 
#' (U.S. Geological Survey, 2012c).
#' Not all statistics are available at every gage.
#' @param sdate is the start date of the time series, in the format 
#' yyyy-mm-dd, optional.
#' @param edate is the end date of the time series, in the format yyyy-mm-dd, 
#' optional.
#' @return a data frame containing daily streamflow or other hydrologic data 
#' for the site specified during the dates specified (inclusive).  The USGS 
#' parameter code, code, and the statistics code, stat, are attributes of the
#' data frame.
#' @references 
#' U.S. Geological Survey, 2012a, National Water Information System: Mapper, 
#' accessed September 7, 2012, at 
#' \url{http://wdr.water.usgs.gov/nwisgmap/}.
#' 
#' U.S. Geological Survey, 2012b, Parameter code definition, 
#' National Water Information System: Web Interface, accessed September 7, 
#' 2012, at \url{http://nwis.waterdata.usgs.gov/usa/nwis/pmcodes}.
#' 
#' U.S. Geological Survey, 2012c, Stat codes (stat_cd), 
#' National Water Information System: Web Interface, accessed September 7, 
#' 2012, at 
#' \url{http://nwis.waterdata.usgs.gov/nwis/help/?read_file=stat&format=table}.
#' 
#' U.S. Geological Survey, 2012d, USGS daily values site web service: 
#' REST Web Services, accessed September 7, 2012, at 
#' \url{http://waterservices.usgs.gov/rest/DV-Service.html}.
#' 
#' U.S. Geological Survey, 2012e, USGS site web service: 
#' REST Web Services, accessed September 7, 2012, at 
#' \url{http://waterservices.usgs.gov/rest/Site-Service.html}.
#' 
#' U.S. Geological Survey, 2012f, USGS surface-water daily data for the Nation: 
#' National Water Information System: Web Interface, accessed September 7, 
#' 2012, at \url{http://waterdata.usgs.gov/nwis/dv/?referred_module=sw}.
#' @keywords ts IO
#' @export
#' @format The returned data frame has the following columns \cr
#' \tabular{lll}{
#' Name \tab Type \tab Description \cr 
#' staid \tab factor \tab USGS station identification number \cr
#' val \tab numeric \tab The value of the hydrologic variable \cr
#' dates \tab Date \tab Date of daily value \cr
#' qualcode \tab factor \tab Qualification code
#' }
#' @examples
#' # import mean daily streamflow for Red River of the North at Fargo, ND
#' q05054000 <- importDVs("05054000", sdate="2000-01-01", edate="2010-12-31")
#' head(q05054000)
#' # additional examples of how to this function follow
#' # import mean daily gage height for Red River of the North at Grand Forks, ND
#' gh05082500 <- importDVs("05082500", code="00065", sdate="2000-01-01", edate="2010-12-31")
#' # import mean daily specific conductance for Red River of the North at Grand Forks, ND
#' sc05082500<- importDVs("05082500", code="00095", sdate="2000-01-01", edate="2010-12-31")
#' # import mean daily water temperature for Red River of the North at Fargo, ND
#' temp05054000<- importDVs("05054000", code="00010", sdate="2000-01-01", edate="2010-12-31")
#' # import median daily pH for Red River of the North at Fargo, ND
#' pH05054000<- importDVs("05054000", code="00400", stat="00008", 
#' sdate="2000-01-01", edate="2010-12-31")
#' # examine the attributes of the data frame to show that the parameter code 
#' # and statistics code are saved with the data frame
#' attributes(pH05054000)[c("code","stat")]
#' # import mean daily oxygen for Red River of the North at Fargo, ND
#' do05054000 <- importDVs("05054000", code="00300", sdate="2000-01-01", edate="2010-12-31")
#' # import mean daily turbidity for Red River of the North at Fargo, ND
#' turb05054000 <- importDVs("05054000", code="63680", sdate="2000-01-01", edate="2010-12-31")
importDVs <- function(staid, code="00060", stat="00003", sdate="1851-01-01", 
                      edate=as.Date(Sys.Date(), format="%Y-%m-%d")) {
  if (is.character(staid) == FALSE ) stop("staid needs to have quotes around it")
  if (nchar(staid) < 8) stop ("staid must be at least 8 characters")
  base_url <- "http://waterservices.usgs.gov/nwis/dv?"
  url <- paste(base_url, "site=", staid, "&parameterCd=", code, "&statCd=", 
               stat, sep = "")
  url <- paste(url, "&startDt=", sdate, "&endDt=", edate, sep="")
  doc <- xmlTreeParse(url, getDTD = FALSE, useInternalNodes=TRUE)
  # Get everything in the main (root) element:
  r <- xmlRoot(doc)
  # Put all the measured values in a vector
  i <- 1
  val <- vector(mode="numeric", length=1)
  while (xmlName(r[[2]][[3]][[i]])=="value") {
    val[i] <- as.numeric(xmlValue(r[[2]][[3]][[i]]))
	  i <- i + 1
  }

  # Put all of the attributes in a Attribute list (there are 2 attributes, 
  # qualifiers and dateTime)
  Attribute <- xmlApply(r[[2]][[3]], xmlAttrs)
 
  # Get the number of data values
  N <- length(val)

  # Get the default NoDataValue
  NoDataValue <- xmlValue(r[["timeSeries"]][["variable"]][["NoDataValue"]])
  NoDataValue <- as.integer(NoDataValue)
  dates <- vector(mode="character", length=1)
  qualcode <- vector(mode="character", length=1)
  if ( N > 1 ) { 
  for (z in 1:N) {
	  dates[z] <- as.character(strsplit(Attribute[z][[1]][[2]], "T")[[1]][1])
    qualcode[z] <- Attribute[z][[1]][[1]]
  }
  dates <- as.Date(dates, "%Y-%m-%d")
  df <- data.frame(staid, val, dates, qualcode)
  beginDate <- df$dates[1]
  endDate <- df$dates[dim(df)[[1]]]
  myDates <- as.data.frame(seq.Date(beginDate, endDate, by=1))
  dimnames(myDates)[[2]][1]<-"dates"
  ndays<-dim(myDates)[1]
  nobs<-dim(df)[1]
  if ( nobs < ndays ) {
    sitedat<-df
    fixedData<-merge(myDates, sitedat, all.x=TRUE)
    fixedData$staid<-sitedat$staid[1]
    fixedData<-fixedData[,c("staid", "val", "dates", "qualcode")]
    df<-fixedData 
  }
  } 
  else { 
    df <- data.frame(staid=character(0), val=numeric(0), dates=character(0), 
                     qualcode=character(0)) 
    my.message<-paste("No data returned for site", staid, "parameter code", code, 
                      "statistics code", stat, sdate, "to", edate, sep=" ")
    message(my.message)
  }
  attributes(df)$code<-code
  attributes(df)$stat<-stat
  df
}

#' Function to plot hydrologic times series.  
#' Will plot more than one site at a time.
#'
#' @name plotParam
#' @title Plot Streamflow and Continous Water-Quality Data
#' @param data is the data frame in the foramt of that returned by 
#' \link{importDVs}.
#' @param metric USGS streamflow data are usually in cubic feet per second;   
#' however it may be converted to cubic meters per second for publication.  
#' Likewise, gage height is usually in feet, but could be converted to 
#' meters.  The metric argument only has an effect on streamflow and gage 
#' height.
#' @param logscale is a logical indicating whether or not the y-scale should be 
#' log 10.  Streamflow generally is plotted with a log scale and this only has 
#' an effect on the plotting of streamflow data.
#' @param ylabel optionally allows user to pass a y-axis label.
#' @param ... further arguments to be passed to plotting method (see \link{par}). 
#' (see \link{xyplot}).
#' @return a lattice plot 
#' @export
#' @examples 
#' data(exampleWaterData)
#' plotParam(misQ05054000, code="00060", stat="00003", logscale=TRUE)
#' plotParam(misQ05054000, code=attributes(misQ05054000)$code, 
#' stat=attributes(misQ05054000)$stat, logscale=TRUE)
#' @keywords hplot ts univar
plotParam<-function(data, logscale=FALSE, metric=FALSE, ylabel=NULL, ...) {
  if (missing(ylabel) ) {
    if ( is.null(attributes(data)$stat) | is.null(attributes(data)$code) ) {
      stop("The data frame needs to have either stat and code attributes or 
           the ylabel argument needs to be used.")
    } else {
      stat<-attributes(data)$stat
      code<-attributes(data)$code
      if (stat=="00003") { stat.txt <- "Daily mean" } else
        if (stat=="00001") { stat.txt <- "Daily maximum" } else
          if (stat=="00002") { stat.txt <- "Daily minimum" } else
            if (stat=="00008") { stat.txt <- "Daily median" } else {
              message("Unknown stat code, a label should be passed using the 
                      ylabel argument.")
            }
      
      if (code=="00060") {
        if ( logscale== TRUE ) {
          if (metric==TRUE) {
            my.ylab <- paste(stat.txt,"streamflow, cubic meters per second", 
                             sep=" ")
            my.plot<-xyplot((val*0.0283)~dates | staid, data=data, typ="l", 
                            scales=list(x = list(tck = -1), y=list(log=TRUE, 
                                                                   tck = -1)), 
                            ylab=my.ylab, xlab="",
                            yscale.components = yscale.components.log10ticks,
                            ...)
          }  
          else if (metric=="FALSE") {
            my.ylab<-paste(stat.txt,"streamflow, cubic feet per second", 
                           sep=" ")
            my.plot<-xyplot(val~dates|staid, data=data, typ="l", 
                            scales=list(x = list(tck = -1), y=list(log=logscale, 
                                                                   tck = -1)), 
                            ylab=my.ylab, xlab="", 
                            yscale.components = yscale.components.log10ticks,
                            ...) 
          } 
          else {
            stop("metric must be TRUE for cubic meters per second or FALSE
             for cubic feet per second")
          }
        }  
        if ( logscale==FALSE ) {
          if (metric==TRUE) {
            my.ylab <- paste(stat.txt,"streamflow, cubic meters per second", 
                             sep=" ")
            my.plot<-my.plot<-xyplot((val*0.0283)~dates | staid, data=data, 
                                     typ="l", scales=list(x = list(tck = -1), 
                                                          y=list(tck = -1)), 
                                     ylab=my.ylab, xlab="", ...)
          }  
          else if (metric=="FALSE") {
            my.ylab<-paste(stat.txt,"streamflow, cubic feet per second", 
                           sep=" ")
            my.plot<-xyplot(val~dates|staid, data=data, typ="l", 
                            scales=list(x = list(tck = -1), y=list(tck = -1)), 
                            ylab=my.ylab, xlab="", ...) 
          } 
          else {
            stop("metric must be TRUE for cubic meters per second or FALSE
           for cubic feet per second")
          }
        }
      }
      if (code=="00065") {
        if (metric==FALSE) {
          my.ylab <- paste(stat.txt, "gage height, feet", sep=" ")
          my.plot<-xyplot(val~dates|staid, data=data, typ="l", 
                          scales=list(x = list(tck = -1), y=list(tck = -1)), 
                          ylab=my.ylab, xlab="", ...)  
        }
        else if (metric==TRUE) {
          my.ylab <- paste(stat.txt, "gage height, meters", sep=" ")
          my.plot<-xyplot((val*0.3048)~dates|staid, data=data, typ="l", 
                          scales=list(x = list(tck = -1), y=list(tck = -1)), 
                          ylab=my.ylab, xlab="", ...)  
        }
        else {
          stop("metric must be TRUE for meters or FALSE for feet.")
        }
      }
      if (code=="00095") {
        my.ylab<-paste(stat.txt, 
                       "specific conductance, water,\nunfiltered, microsiemens per centimeter\nat 25 degrees Celsius", 
                       sep=" ")
        my.plot<-xyplot(val~dates|staid, data=data, typ="l", 
                        scales=list(x = list(tck = -1), y=list(tck = -1)), 
                        ylab=my.ylab, xlab="", ...) 
      }
      if (code=="00010") {
        my.ylab<-paste(stat.txt, "temperature, water, degrees Celsius", sep=" ")
        my.plot<-xyplot(val~dates|staid, data=data, typ="l", 
                        scales=list(x = list(tck = -1), y=list(tck = -1)), 
                        ylab=my.ylab, xlab="", ...)
      }
      if (code=="00300") {
        my.ylab <- paste(stat.txt, 
                         "dissolved oxygen, water, unfiltered, milligrams per liter", 
                         sep=" ")
        my.plot<-xyplot(val~dates|staid, data=data, typ="l", 
                        scales=list(x = list(tck = -1), y=list(tck = -1)), 
                        ylab=my.ylab, xlab="", ...)  
      }
      if (code=="00400") {
        my.ylab <- paste(stat.txt, 
                         "pH, water, unfiltered, field, standard units", 
                         sep=" ")
        my.plot<-xyplot(val~dates|staid, data=data, typ="l", 
                        scales=list(x = list(tck = -1), y=list(tck = -1)), 
                        ylab=my.ylab, xlab="", ...)  
      }
      if (code=="63680") {
        my.ylab <- paste(stat.txt, 
                         "turbidity, water, unfiltered,\nmonochrome near infra-red LED light,\n780-900 nm, detection angle 90 +/ -2.5 degrees,\nformazin nephelometric units (FNU)", 
                         sep=" ")
        my.plot<-xyplot(val~dates|staid, data=data, typ="l", 
                        scales=list(x = list(tck = -1), y=list(tck = -1)), 
                        ylab=my.ylab, xlab="", ...) 
      }
    }  
  } else {
    my.ylab <- ylabel
    if ( logscale==TRUE ) {
      my.plot<-xyplot(val~dates | staid, data=data, typ="l", 
                      scales=list(x = list(tck = -1), y=list(log=TRUE, 
                                                               tck = -1)), 
                      ylab=my.ylab, xlab="", 
                      yscale.components = yscale.components.log10ticks, ...)
    } else {
      my.plot<-xyplot(val~dates | staid, data=data, typ="l", 
                      scales=list(x = list(tck = -1), y=list(log=FALSE, 
                                                             tck = -1)), 
                      ylab=my.ylab, xlab="")
    }   
  }
  my.plot
}

#' Function that returns USGS Daily Values Site Service URL for troubleshooting or 
#' building a URL for other purposes.
#'
#' @name tellMeURL
#' @title USGS Daily Values Site Service URL
#' @param staid is the USGS site identification number, which 
#' is usually eight digits long, but can be longer.  Users may search for 
#' surface-water sites and obtain station identification numbers using the 
#' USGS Site Web Service, 
#' \url{http://waterservices.usgs.gov/rest/Site-Service.html} (U.S. Geological 
#' Survey, 2012d); using the National Water Information System: Mapper, 
#' \url{http://wdr.water.usgs.gov/nwisgmap/} (U.S. Geological Survey, 2012a); 
#' or using the National Water Information System: Web Interface to daily 
#' surface-water data, 
#' \url{http://waterdata.usgs.gov/nwis/dv/?referred_module=sw} (U.S. 
#' Geological Survey, 2012e).  The site identification number needs to be 
#' entered as a character, that is in quotes, because many USGS streamgage 
#' numbers begin with zero and the leading zero is necessary.
#' @param code is the USGS parameter code, a 5-digit number 
#' used in the USGS computerized data system, National Water 
#' Information System (NWIS), to uniquely identify a specific hydrologic 
#' property or constituent.  A list of paramater codes is available at 
#' \url{http://nwis.waterdata.usgs.gov/usa/nwis/pmcodes} (U.S. Geological 
#' Survey, 2012b).
#' @param stat is the USGS statistics code, a 5-digit number 
#' used in the USGS computerized data system, NWIS, to uniquely identify 
#' specific statistics, such as daily mean, daily maximum, and daily minimum.  
#' The default,  00003,  is the mean daily value.  A list of statistics codes 
#' is available at 
#' \url{http://nwis.waterdata.usgs.gov/nwis/help/?read_file=stat&format=table} 
#' (U.S. Geological Survey, 2012c).
#' Not all statistics are available at every gage.
#' @param sdate is the start date of the time series, in the format yyyy-mm-dd.
#' @param edate is the end date of the time series, in the format yyyy-mm-dd.
#' @keywords utilities
#' @export
#' @references
#' U.S. Geological Survey, 2012a, National Water Information System: Mapper, 
#' accessed September 7, 2012, at 
#' \url{http://wdr.water.usgs.gov/nwisgmap/}.
#' 
#' U.S. Geological Survey, 2012b, Parameter code definition, 
#' National Water Information System: Web Interface, accessed September 7, 
#' 2012, at \url{http://nwis.waterdata.usgs.gov/usa/nwis/pmcodes}.
#' 
#' U.S. Geological Survey, 2012c, Stat codes (stat_cd), 
#' National Water Information System: Web Interface, accessed September 7, 
#' 2012, at 
#' \url{http://nwis.waterdata.usgs.gov/nwis/help/?read_file=stat&format=table}.
#' 
#' U.S. Geological Survey, 2012d, USGS site web service: 
#' REST Web Services, accessed September 7, 2012, at 
#' \url{http://waterservices.usgs.gov/rest/Site-Service.html}.
#' 
#' U.S. Geological Survey, 2012e, USGS surface-water daily data for the Nation: 
#' National Water Information System: Web Interface, accessed September 7, 
#' 2012, at \url{http://waterdata.usgs.gov/nwis/dv/?referred_module=sw}.
#' @return URL for USGS data
#' @examples 
#' tellMeURL("05054000", code="00060", stat="00003", sdate="2000-01-01", 
#'  edate=as.Date(Sys.Date(), format="%Y-%m-%d"))
tellMeURL <- function(staid, code="00060", stat="00003", sdate="1851-01-01",  
                      edate=as.Date(Sys.Date(), format="%Y-%m-%d")) {
  if (is.character(staid) == FALSE ) stop("staid needs to have quotes around it")
  if (nchar(staid) < 8) stop ("staid must be at least 8 characters")
  base_url <- "http://waterservices.usgs.gov/nwis/dv?"
  url <- paste(base_url, "site=", staid, "&parameterCd=", code, "&statCd=", 
               stat, sep = "")
  url <- paste(url, "&startDt=", sdate, sep="")
  url <- paste(url, "&endDt=", edate, sep="")
  url
}

#' Function to retrieve information about a USGS streamgage site
#' 
#' This provides some limited metadata about the USGS streamgage site.
#' @name siteInfo
#' @title Retrieve site information
#' @param staid is a character vector containing USGS site
#' identification number(s).  USGS site numbers are usually eight digits long, 
#' but can be longer.  Users may search for surface-water sites and obtain 
#' station identification numbers using the USGS Site Web Service, 
#' \url{http://waterservices.usgs.gov/rest/Site-Service.html} (U.S. Geological 
#' Survey, 2012b); using the National Water Information System: Mapper, 
#' \url{http://wdr.water.usgs.gov/nwisgmap/} (U.S. Geological Survey, 2012a); 
#' or using the National Water Information System: Web Interface to daily 
#' surface-water data, 
#' \url{http://waterdata.usgs.gov/nwis/dv/?referred_module=sw} (U.S. 
#' Geological Survey, 2012c).  The site identification number needs to be 
#' entered as a character, that is in quotes, because many USGS streamgage 
#' numbers begin with zero and the leading zero is necessary.
#' @keywords datagen
#' @references
#' U.S. Geological Survey, 2012a, National Water Information System: Mapper, 
#' accessed September 7, 2012, at 
#' \url{http://wdr.water.usgs.gov/nwisgmap/}.
#' 
#' U.S. Geological Survey, 2012b, USGS site web service: 
#' REST Web Services, accessed September 7, 2012, at 
#' \url{http://waterservices.usgs.gov/rest/Site-Service.html}.
#' 
#' U.S. Geological Survey, 2012c, USGS surface-water daily data for the Nation: 
#' National Water Information System: Web Interface, accessed September 7, 
#' 2012, at \url{http://waterdata.usgs.gov/nwis/dv/?referred_module=sw}.
#' @export
#' @return a data frame containing the station identification number(s), the 
#' USGS streamgage name(s), the decimal latitude(s), and decimal longitude(s).
#' @format a data frame with the following columns:
#'   \tabular{llll}{
#' Name \tab Type \tab Description \cr 
#' staid \tab factor \tab USGS station identification number \cr
#' staname \tab character \tab USGS station name \cr
#' lat \tab numeric \tab Decimal latitude \cr
#' lng \tab numeric \tab Decimal longitude
#' }
#' @note Information retrieved using this function can be used to create a map of
#' multiple streamgage sites---see package vignette.
#' @examples
#' staInfo <- siteInfo("05054000")
#' staInfo
#' staInfo <- siteInfo(c("05054000", "05082500", "06342500"))
#' staInfo
#' # a list with an invalid station identification number
#' staInfo <- siteInfo(c("05054000", "05082500", "06342501"))
siteInfo<-function(staid) {
  staname <- vector(mode="character", length=0)
  lat <- vector(mode="numeric", length=0)
  lng <- vector(mode="numeric", length=0)
  for ( i in 1:length(staid) ) {
    if (is.character(staid[i]) == FALSE ) stop("staid needs to have quotes 
                                               around it")
    if (nchar(staid[i]) < 8) stop ("staid must be at least 8 characters")
    base_url <-"http://waterservices.usgs.gov/nwis/site?format=mapper&sites="
    url <- paste(base_url,staid[i],
                 "&siteOutput=expanded&seriesCatalogOutput=true&outputDataTypeCd=all", 
                 sep = "")
    my.er<- FALSE
    doc <- try(xmlTreeParse(url, getDTD = FALSE, useInternalNodes=TRUE), 
               silent=TRUE)
    if ( class(doc)[1]=="try-error") {
      message("Unable to retrieve XML for site ", staid[i], ". Check site id.")
      staname[i] <- "Unable to retrieve"
      lat[i] <- ""
      lng[i] <- ""
      next
    }
    else {    
      # Get everything in the main ('root') element:
      r <- xmlRoot(doc)
      r.attrs <- xmlApply(r[[1]], xmlAttrs)
      staname[i] <- r.attrs$site["sna"]
      lat[i] <- r.attrs$site["lat"]
      lng[i] <- r.attrs$site["lng"]
    }
  }
  df <- data.frame(cbind(staid=staid, staname=staname, lat=lat, lng=lng), 
                   stringsAsFactors=FALSE)
  df$staname <- as.character(df$staname)
  df$lat <- as.numeric(lat)
  df$lng <- as.numeric(lng)
  df
}

#' Function that returns USGS Site Information Service URL for troubleshooting or 
#' building a URL for other purposes.
#'
#' @name tellMeSiteURL
#' @title USGS Site Information Service URL
#' @param staid is the USGS site identification number, which 
#' is usually eight digits long, but can be longer.  Users may search for 
#' surface-water sites and obtain station identification numbers using the 
#' USGS Site Web Service, 
#' \url{http://waterservices.usgs.gov/rest/Site-Service.html} (U.S. Geological 
#' Survey, 2012b); using the National Water Information System: Mapper, 
#' \url{http://wdr.water.usgs.gov/nwisgmap/} (U.S. Geological Survey, 2012a); 
#' or using the National Water Information System: Web Interface to daily 
#' surface-water data, 
#' \url{http://waterdata.usgs.gov/nwis/dv/?referred_module=sw} (U.S. 
#' Geological Survey, 2012c).  The site identification number needs to be 
#' entered as a character, that is in quotes, because many USGS streamgage 
#' numbers begin with zero and the leading zero is necessary.
#' @keywords utilities
#' @return URL for USGS site information
#' @export
#' @references
#' U.S. Geological Survey, 2012a, National Water Information System: Mapper, 
#' accessed September 7, 2012, at 
#' \url{http://wdr.water.usgs.gov/nwisgmap/}.
#' 
#' U.S. Geological Survey, 2012b, USGS site web service: 
#' REST Web Services, accessed September 7, 2012, at 
#' \url{http://waterservices.usgs.gov/rest/Site-Service.html}.
#' 
#' U.S. Geological Survey, 2012c, USGS surface-water daily data for the Nation: 
#' National Water Information System: Web Interface, accessed September 7, 
#' 2012, at \url{http://waterdata.usgs.gov/nwis/dv/?referred_module=sw}.
#' @examples 
#' tellMeSiteURL("05054000")
tellMeSiteURL <- function(staid) {
if (is.character(staid) == FALSE ) stop("Station number needs to have quotes around it")
  if (nchar(staid) < 8) stop ("Station number must be at least 8 characters")
  base_url <-"http://waterservices.usgs.gov/nwis/site?format=mapper&sites="
  url <- paste(base_url,staid,
               "&siteOutput=expanded&seriesCatalogOutput=true&outputDataTypeCd=all", 
               sep = "")
  url
}
