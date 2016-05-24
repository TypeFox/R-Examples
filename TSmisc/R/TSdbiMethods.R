########################################################################
##                 histQuote
########################################################################

histQuote <- function(...) {
  drv <- "histQuote"
  attr(drv, "package") <- "TShistQuote"
  new("histQuoteDriver", Id = drv)
  }

####### some kludges to make this look like DBI. ######
# for this require("DBI")

setClass("histQuoteDriver", contains=c("DBIDriver"), slots=c(Id = "character")) 

setClass("histQuoteConnection", contains=c("DBIConnection", "histQuoteDriver"),
   slots=c(dbname="character") )

setMethod("dbConnect", signature(drv="histQuoteDriver"), 
     definition=function(drv, dbname, ...) new("histQuoteConnection", drv, dbname=dbname))

# this does nothing but prevent errors if it is called. 
setMethod("dbDisconnect", signature(conn="histQuoteConnection"), 
     definition=function(conn, ...) TRUE)
#######     end kludges   ######

setClass("TShistQuoteConnection", contains=c("histQuoteConnection", "conType","TSdb"))

setMethod("TSconnect",   signature(q="histQuoteConnection", dbname="missing"), 
  definition= function(q, dbname, ...){
   dbname <- q@dbname
   if (is.null(dbname)) stop("dbname must be specified")
   if (dbname == "yahoo") {
      con <- try(url("http://quote.yahoo.com"), silent = TRUE)
      if(inherits(con, "try-error")) 
         stop("Could not establish TShistQuoteConnection to ",  dbname)
      close(con)
      }
   else if (dbname == "oanda") {
      con <- try(url("http://www.oanda.com"),   silent = TRUE)
      if(inherits(con, "try-error")) 
         stop("Could not establish TShistQuoteConnection to ",  dbname)
      close(con)
      }
   else 
      warning(dbname, "not recognized. Connection assumed working, but not tested.")
   
   new("TShistQuoteConnection", q="histQuote", dbname=dbname, hasVintages=FALSE, hasPanels=FALSE) 
   } )


setMethod("TSdates",
  signature(serIDs="character", con="TShistQuoteConnection", vintage="ANY", panel="ANY"),
   definition= function(serIDs, con, vintage=NULL, panel=NULL, ... )  
{  # Indicate  dates for which data is available.
   # This requires retrieving series individually so they are not truncated.
   r <- av <- st <- en <- tb <- NULL
   for (i in 1:length(serIDs))
     {r <- try(TSget(serIDs[i], con), silent = TRUE)

      if(inherits(r, "try-error") ) {
        av <- c(av, FALSE)
	st <- append(st, list(NA))
	en <- append(en, list(NA))
	tb <- rbind(tb, NA)
	}
      else  {
        av <- c(av, TRUE)
        st <- append(st, list(tfstart(r)))
        en <- append(en, list(tfend(r)))
        tb <- rbind(tb,tffrequency(r))
        }
      }
  r <- serIDs
  attr(r, "TSdates") <- av
  attr(r, "start") <- st
  attr(r, "end")   <- en
  attr(r, "frequency")   <- tb
  class(r) <- "TSdates"
  r
} )


setMethod("TSget",     signature(serIDs="character", con="TShistQuoteConnection"),
   definition= function(serIDs, con, TSrepresentation=options()$TSrepresentation,
       tf=NULL, start=tfstart(tf), end=tfend(tf),
       names=serIDs, quote = "Close", quiet=TRUE, repeat.try=3, ...){ 

    mat <- desc <- NULL
    # recycle serIDs and quote to matching lengths
    # argument 'quote' ignored for provider 'oanda'
    if (con@dbname == "yahoo") {
       if (length(quote) < length(serIDs))
           quote  <- rep(quote,  length.out=length(serIDs))
       if (length(quote) > length(serIDs))
           serIDs <- rep(serIDs, length.out=length(quote))
       }
    # this is ugly because missing or null args cannot be passed
    args <- list( quiet=quiet, provider = con@dbname, retclass = "zoo")
    args <- if (is.null(start) & is.null(end)) append(args, list(...))
            else if (is.null(start)  ) append(args, list(end=end, ...))
            else if (is.null(end)  )   append(args, list(start=start, ...))
            else                       append(args, list(start=start, end=end, ...) )
    for (i in seq(length(serIDs))) {
       argsi <- if (con@dbname == "yahoo")       
	        append(list(instrument=serIDs[i], quote=quote[i]), args)
              else if (con@dbname == "oanda") 
	        append(list(instrument=serIDs[i]),  args)
       for (rpt in seq(repeat.try)) {
           r <- try(do.call("get.hist.quote", argsi))
	   if (!inherits(r , "try-error")) break
	   }
       if (inherits(r , "try-error")) stop(r)
       if (is.character(r)) stop(r)
       mat <- tbind(mat, r)
       desc <- c(desc, paste(serIDs[i], quote[i], collapse=" "))
       }
    if (NCOL(mat) != length(serIDs)) stop("Error retrieving series", serIDs) 

    mat <- tfwindow(mat, tf=tf, start=start, end=end)
  
    mat <- tframePlus::changeTSrepresentation(mat, TSrepresentation)

    if(all(quote == quote[1])) TSrefperiod(mat) <- quote[1]

    seriesNames(mat) <- names

    TSmeta(mat) <- new("TSmeta", serIDs=serIDs,  dbname=con@dbname, 
        hasVintages=con@hasVintages, hasPanels=con@hasPanels,
  	conType=class(con), DateStamp= Sys.time(), 
	TSdoc=paste(desc, " from ", con@dbname, "retrieved ", Sys.time()),
	TSdescription=paste(desc, " from ", con@dbname),
	TSlabel=desc, 
	TSsource= con@dbname)
    mat
    } )


#setMethod("TSput",     signature(x="ANY", serIDs="character", con="TShistQuoteConnection"),
#   definition= function(x, serIDs=seriesNames(data), con, ...)   
#    "TSput for TShistQuote connection not supported." )

setMethod("TSdescription",   signature(x="character", con="TShistQuoteConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSdescription for TShistQuote connection not supported." )


setMethod("TSdoc",   signature(x="character", con="TShistQuoteConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSdoc for TShistQuote connection not supported." )

setMethod("TSlabel",   signature(x="character", con="TShistQuoteConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSlabel for TShistQuote connection not supported." )


setMethod("TSsource",   signature(x="character", con="TShistQuoteConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        con@dbname)


########################################################################
##                 getSymbol
########################################################################

getSymbol <- function(...) {
  drv <- "getSymbol"
  attr(drv, "package") <- "TSgetSymbol"
  new("getSymbolDriver", Id = drv)
  }

####### some kludges to make this look like DBI. ######
# for these require("DBI")

setClass("getSymbolDriver", contains="DBIDriver", slots=c(Id = "character")) 

setClass("getSymbolConnection", contains=c("DBIConnection", "getSymbolDriver"),
   slots=c(dbname="character") )

setMethod("dbConnect", signature(drv="getSymbolDriver"), 
     definition=function(drv, dbname, ...) new("getSymbolConnection", drv, dbname=dbname))

# this does nothing but prevent errors if it is called. 
setMethod("dbDisconnect", signature(conn="getSymbolConnection"), 
     definition=function(conn,...) TRUE)

#######     end kludges   ######

setClass("TSgetSymbolConnection", contains=c("getSymbolConnection","conType", "TSdb")) 
 
setMethod("TSconnect",   signature(q="getSymbolConnection", dbname="missing"),
  definition= function(q, dbname, ...){
   dbname <- q@dbname
   if (dbname == "FRED") {
      #there could be a better test
      con <- try(quantmod::getSymbols('CPIAUCNS',src='FRED'), silent = TRUE)
      if(inherits(con, "try-error")) 
         stop("Could not establish TSgetSymbolConnection to ",  dbname)
      #close(con)
      }
   else if (dbname == "yahoo") {
      #this breaks if the symbol disappears, so it is more trouble than value
      # a better test would be good
      #con <- try(quantmod::getSymbols('QQQQ',src='yahoo'), silent = TRUE)
      #if(inherits(con, "try-error")) 
      #   stop("Could not establish TSgetSymbolConnection to ",  dbname)
      ##close(con)
      }
   else 
      warning(dbname, "not recognized. Connection assumed working, but not tested.")
   
   new("TSgetSymbolConnection", dbname=dbname,
          hasVintages=FALSE, hasPanels=FALSE) 
   } )


setMethod("TSdates",
  signature(serIDs="character", con="TSgetSymbolConnection", vintage="ANY", panel="ANY"),
   definition= function(serIDs, con, vintage=NULL, panel=NULL, ... )  
{  # Indicate  dates for which data is available.
   # This requires retrieving series individually so they are not truncated.
   r <- av <- st <- en <- tb <- NULL
   for (i in 1:length(serIDs))
     {r <- try(TSget(serIDs[i], con), silent = TRUE)

      if(inherits(r, "try-error") ) {
        av <- c(av, FALSE)
	st <- append(st, list(NA))
	en <- append(en, list(NA))
	tb <- rbind(tb, NA)
	}
      else  {
        av <- c(av, TRUE)
        st <- append(st, list(tfstart(r)))
        en <- append(en, list(tfend(r)))
        tb <- rbind(tb,tffrequency(r))
        }
      }
  r <- serIDs
  attr(r, "TSdates") <- av
  attr(r, "start") <- st
  attr(r, "end")   <- en
  attr(r, "frequency")   <- tb
  class(r) <- "TSdates"
  r
} )

#trace("TSget", browser, exit=browser, signature = c(serIDs="character", #con="TSgetSymbolConnection"))

setMethod("TSget",     signature(serIDs="character", con="TSgetSymbolConnection"),
   definition= function(serIDs, con, TSrepresentation=options()$TSrepresentation,
       tf=NULL, start=tfstart(tf), end=tfend(tf),
       names=serIDs, quote = NULL, 
       quiet=TRUE, repeat.try=3, ...){ 

    #This TSget passes of TSrepresentation to getSymbols rather than 
    #  using tframePlus::changeTSrepresentation
    if (is.null(TSrepresentation)) {
       default <- TRUE
       TSrepresentation <- "zoo"
       }
    else default <- FALSE
    
    if (! TSrepresentation %in% c("ts", "its", "zoo", "xts", "timeSeries"))
       stop(TSrepresentation, " time series class not supported.")

    mat <- desc <- NULL
    # recycle serIDs and quote to matching lengths
    # argument 'quote' ignored for provider 'oanda'
    # if quote is null then HLOC will be retained
    if (con@dbname == "yahoo" && !is.null(quote)) {
        if (length(quote) < length(serIDs))
            quote  <- rep(quote,  length.out=length(serIDs))
        if (length(quote) > length(serIDs))
            serIDs <- rep(serIDs, length.out=length(quote))
        }
    
    args <- list(src = con@dbname, return.class=TSrepresentation,
                 auto.assign=FALSE)
    
    if (con@dbname == "yahoo" )
       args <- if (is.null(start) & is.null(end)) append(args, list(...))
            else if (is.null(start)  ) append(args, list(to=end, ...))
            else if (is.null(end)  )   append(args, list(from=start, ...))
            else         append(args, list(from=start, to=end, ...) )

    for (i in seq(length(serIDs))) {
       argsi <- append(list(serIDs[i]),  args)
       for (rpt in seq(repeat.try)) {
           # quantmod::getSymbols
           r <- try(do.call("getSymbols", argsi), silent=quiet)
	   if (!inherits(r , "try-error")) break
	   }
       if (inherits(r , "try-error")) stop("series not retrieved:", r)
       if (is.character(r)) stop("series not retrieved:", r)
       #TSrefperiod(r) <- quote[i]
       if (!is.null(quote)){
          id <- toupper(sub("^", "", serIDs[i], fixed=TRUE))
	  r <- r[, paste(id,".", quote[i], sep="")]
	  }
       mat <- tbind(mat, r)
       desc <- c(desc, paste(serIDs[i], quote[i], collapse=" "))
       }
    #if (NCOL(mat) != length(serIDs)) stop("Error retrieving series", serIDs)
    #  yahoo connections return high, low , ... 
    if (NCOL(mat) != length(serIDs)) names <- seriesNames(mat) 
    
    st <- as.POSIXlt(start(mat)) #POSIXlt as return for zoo
    if (default) {
        if(xts::periodicity(mat)$scale == "monthly")
	   mat <- ts(mat, frequency=12,start=c(1900+st$year, 1+st$mon))
        else if(xts::periodicity(mat)$scale == "quarterly")
	   mat <- ts(mat, frequency=4, start=c(1900+st$year, 1+(st$mon)/3))
        else if(xts::periodicity(mat)$scale == "yearly")  
	   mat <- ts(mat, frequency=1, start=c(1900+st$year, 1))
	}

    if (con@dbname != "yahoo" )
        mat <- tfwindow(mat, tf=tf, start=start, end=end)

    seriesNames(mat) <- names

    TSmeta(mat) <- new("TSmeta", serIDs=serIDs,  dbname=con@dbname, 
        hasVintages=con@hasVintages, hasPanels=con@hasPanels,
  	conType=class(con), DateStamp= Sys.time(), 
	TSdoc=paste(desc, " from ", con@dbname, "retrieved ", Sys.time()),
	TSdescription=paste(desc, " from ", con@dbname),
	TSlabel=desc, 
	TSsource= (if("yahoo" == con@dbname) "yahoo" 
	      else if("FRED" == con@dbname) "Federal Reserve Bank of St. Louis"
	      else con@dbname )
	) 
    mat
    } 
    )


#setMethod("TSput",     signature(x="ANY", serIDs="character", con="TSgetSymbolConnection"),
#   definition= function(x, serIDs=seriesNames(data), con, ...)   
#    "TSput for TSgetSymbol connection not supported." )

setMethod("TSdescription",   signature(x="character", con="TSgetSymbolConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSdescription for TSgetSymbol connection not supported." )


setMethod("TSdoc",   signature(x="character", con="TSgetSymbolConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSdoc for TSgetSymbol connection not supported." )

setMethod("TSlabel",   signature(x="character", con="TSgetSymbolConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSlabel for TSgetSymbol connection not supported." )

setMethod("TSsource",   signature(x="character", con="TSgetSymbolConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        if("yahoo" == con@dbname) "yahoo" 
	      else if("FRED" == con@dbname) "Federal Reserve Bank of St. Louis"
	      else "unspecified"  )

#######################################################################
##                 xls
########################################################################
xls <- function(...) {
  drv <- "xls"
  attr(drv, "package") <- "TSxls"
  new("xlsDriver", Id = drv)
  }

####### some kludges to make this look like DBI. ######
# for this require("DBI")

setClass("xlsDriver", contains=c("DBIDriver"), slots=c(Id = "character")) 

setClass("xlsConnection", contains=c("DBIConnection", "xlsDriver"),
   slots=c(dbname="character") )

setMethod("dbConnect", signature(drv="xlsDriver"), 
     definition=function(drv, dbname, ...) 
                   new("xlsConnection", drv, dbname=dbname))

# this does nothing but prevent errors if it is called. 
setMethod("dbDisconnect", signature(conn="xlsConnection"), 
     definition=function(conn,...) TRUE)
#######     end kludges   ######


setClass("TSxlsConnection", contains=c("xlsConnection", "conType","TSdb"),
   slots= c(url="character", data="matrix", ids="character", 
        dates="character", names="character", description="character",
	source="character", tsrepresentation = "function") 
   )

# ... is passed  by default TSconnect (in TSdbi) to both the
#  the driver  and also to  TSconnect con signature methods.
# Here map is used by TSconnect but in most other cases it is used by
#  the driver, so this would be more similar to other packages if it
#  was stored by the driver (and then would not need to be passed to both).

setMethod("TSconnect",   signature(q="xlsConnection", dbname="missing"),
  definition= function(q, dbname, 
     map=list(ids, data, dates, names=NULL, description=NULL, sheet=1,
              tsrepresentation = function(data,dates){
		       zoo(data, as.Date(dates))}), ...){
   #  user / password / host  for future consideration
   dbname <- q@dbname 

   sheet <- if (is.null(map$sheet)) 1 else map$sheet

   if (file.exists(dbname)) {
      file <- dbname
      url <- ""
      }
   else{
     url <- dbname
     file <- tempfile()
     on.exit(unlink(file) )
     zz <- try(download.file(url, file, quiet = FALSE, mode = "wb",
                   cacheOK = TRUE),  silent=TRUE)
     #or url(url)

     if(inherits(zz, "try-error") || (0 != zz)) 
       stop("download.file error, possibly could not find url ",  url,
            " or file ", file)
     }

    if (requireNamespace("gdata", quietly = TRUE)) {
      zz <- try(gdata::read.xls(file, sheet=sheet, blank.lines.skip=FALSE,
           verbose=FALSE), silent=TRUE) #method=c("csv","tsv","tab"), perl="perl")
      if(inherits(zz, "try-error")) 
           stop("Could not read spreedsheet ",  dbname, zz)
      } else warning("gdata is needed for TSxls connection.")

   #NB The first line provides data frame names, so rows are shifted. 
   #   This fixes so matrix corresponds to spreadsheet cells
   z <- rbind(names(zz), as.matrix(zz))

   # translate cell letter range to number indices
   jmap <- function(cols){ 
	st <- unlist(strsplit(sub(":[A-Z]*","",cols),""))
	en <- unlist(strsplit(sub("[A-Z]*:","",cols),""))
	sum(  charmatch(st, LETTERS) * 26^c(0:(length(st)-1))):
	  sum(charmatch(en, LETTERS) * 26^c(0:(length(en)-1)))
	}

   ids   <- z[map$ids$i,  jmap(map$ids$j)] 
   data  <- z[map$data$i, jmap(map$data$j), drop=FALSE]   
   dates <- z[map$dates$i,jmap(map$dates$j)]   
   nm    <- if(is.null(map$names)) NULL else combineRows(
            z[map$names$i,jmap(map$names$j), drop=FALSE]) 
   desc  <- if(is.null(map$description)) NULL else combineRows(
            z[map$description$i,jmap(map$description$j)]) 
   
   #seriesInColumns=TRUE, assuming this for now
   
   z <- dim(data)
   data <- try(as.numeric(data),  silent=TRUE)
   if(inherits(data, "try-error")) 
         stop("Error converting  data to numeric.", data)
   
   data <- array(data, z)

   if(length(dates) != NROW(data))
       stop("length of dates not equal length of series.")
   
   if(length(ids)   != NCOL(data))
       stop("number of ids not equal number of series.")
   
   if(!is.null(nm)) if(length(nm)   != NCOL(data))
       stop("number of names not equal number of series.")

   if(!is.null(desc)) if(length(desc)   != NCOL(data))
       stop("number of descriptions not equal number of series.")

   if(is.null(nm))   nm   <- rep("",NCOL(data))
   if(is.null(desc)) desc <- rep("",NCOL(data))
   
   seriesNames(data ) <- ids
   names(nm)   <- ids
   names(desc) <- ids

   #Adjustments <- c(rep("nsa", 10),rep("sa", 4),rep("nsa", 2)) 
   #Units    <- z[1, 1]  ; names(Units) <- NULL
   #Notes    <- z[2, 1]  ; names(Notes) <- NULL
   #Updated <- z[8, -1]  ; names(Updated) <- NULL# date format error?
   #Source   <- z[9, -1] ; names(Source) <- NULL
 
   # check that tsrepresentation works
   z <- try(map$tsrepresentation(data[,1], dates),  silent=TRUE)
   if(inherits(z, "try-error")) 
         stop("Could not convert data to series using tsrepresentation.",z)
  
   # cache data, etc in con
   # use ids to extract from cache, but give names
   new("TSxlsConnection", dbname=dbname, 
        hasVintages=FALSE, hasPanels=FALSE, url=url,
	data=data,ids=ids,dates=dates, names=nm, description=desc,
	source=dbname,  #this could be better
	tsrepresentation=map$tsrepresentation) 
   } 
   )

setMethod("TSdates",
  signature(serIDs="character", con="TSxlsConnection", vintage="ANY", panel="ANY"),
   definition= function(serIDs, con, vintage=NULL, panel=NULL, ... )  
{  # Indicate  dates for which data is available.
   # This requires retrieving series individually so they are not truncated.
   r <- av <- st <- en <- tb <- NULL
   for (i in 1:length(serIDs))
     {r <- try(TSget(serIDs[i], con), silent = TRUE)

      if(inherits(r, "try-error") ) {
        av <- c(av, FALSE)
	st <- append(st, list(NA))
	en <- append(en, list(NA))
	tb <- rbind(tb, NA)
	}
      else  {
        av <- c(av, TRUE)
        st <- append(st, list(tfstart(r)))
        en <- append(en, list(tfend(r)))
        tb <- rbind(tb,tffrequency(r))
        }
      }
  r <- serIDs
  attr(r, "TSdates") <- av
  attr(r, "start") <- st
  attr(r, "end")   <- en
  attr(r, "frequency")   <- tb
  class(r) <- "TSdates"
  r
} )

setMethod("TSget",     signature(serIDs="character", con="TSxlsConnection"),
   definition=function(serIDs, con, TSrepresentation=options()$TSrepresentation,
       tf=NULL, start=tfstart(tf), end=tfend(tf),
       names=serIDs, ...){ 
    
    # data, ids and dates are cached in con
    mat <- try(con@tsrepresentation(con@data[,serIDs], con@dates),
               silent=TRUE)
    if(inherits(mat, "try-error")) 
         stop("Could not convert data to series using tsrepresentation.",mat)
    # give names rather than id mnemonic 
    seriesNames(mat) <- con@names[serIDs]
    desc <- con@description[serIDs]

    if (NCOL(mat) != length(serIDs)) stop("Error retrieving series", serIDs) 

    mat <- tfwindow(mat, tf=tf, start=start, end=end)
  
    mat <- tframePlus::changeTSrepresentation(mat, TSrepresentation)

    seriesNames(mat) <- names

    TSmeta(mat) <- new("TSmeta", serIDs=serIDs,  dbname=con@dbname, 
        hasVintages=con@hasVintages, hasPanels=con@hasPanels,
  	conType=class(con), DateStamp= Sys.time(), 
	TSdoc=paste(desc, " from ", con@dbname, "retrieved ", Sys.time()),
	TSdescription=paste(desc, " from ", con@dbname),
	TSlabel=desc,
	TSsource=con@dbname ) 
    mat
    } 
    )


#setMethod("TSput",     signature(x="ANY", serIDs="character", con="TSxlsConnection"),
#   definition= function(x, serIDs=seriesNames(data), con, ...)   
#    "TSput for TSxls connection not supported." )

setMethod("TSdescription",   signature(x="character", con="TSxlsConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSdescription for TSxls connection not supported." )


setMethod("TSdoc",   signature(x="character", con="TSxlsConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSdoc for TSxls connection not supported." )

setMethod("TSlabel",   signature(x="character", con="TSxlsConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSlabel for TSxls connection not supported." )

setMethod("TSsource",   signature(x="character", con="TSxlsConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSsource for TSxls connection not supported." )

##### utilities #####

trimAllNA <- function(x, startNAs= TRUE, endNAs= TRUE) UseMethod("trimAllNA") 

trimAllNA.default <- function(x, startNAs= TRUE, endNAs= TRUE)
{# trim NAs from the ends of a ts matrix or vector.
 # (Observations are dropped if all in a given period contain NA.)
 # if startNAs=F then beginning NAs are not trimmed.
 # If endNAs=F   then ending NAs are not trimmed.
 sample <- ! if (is.matrix(x)) apply(is.na(x),1, all) else is.na(x)
 if (!any(sample)) warning("data is empty after triming NAs.")
 s <- if (startNAs) min(time(x)[sample]) else tfstart(x)
 e <- if (endNAs)   max(time(x)[sample]) else tfend(x)
 tfwindow(x, start=s, end=e, warn=FALSE)
}

combineRows <- function(x, i, j, setEmpty=NULL){
  # combine rows of text when it extends over more than a line.
  x[setEmpty] <- ""
  x <- x[i,j, drop=FALSE]
  r <- NULL
  for (ii in seq(NROW(x))) r <- paste(r, x[ii,, drop=FALSE])
  r
  }

#######################################################################
##                 zip
########################################################################
zip <- function(...) {
  drv <- "zip"
  attr(drv, "package") <- "TSzip"
  new("zipDriver", Id = drv)
  }

####### some kludges to make this look like DBI. ######
# for this require("DBI")
setClass("zipDriver", contains=("DBIDriver"), slots=c(Id = "character")) 

setClass("zipConnection", contains=c("DBIConnection", "zipDriver"),
   slots=c(dbname="character") )

setMethod("dbConnect", signature(drv="zipDriver"), 
     definition=function(drv, dbname, ...) 
                   new("zipConnection", drv, dbname=dbname))

# this does nothing but prevent errors if it is called. 
setMethod("dbDisconnect", signature(conn="zipConnection"), 
     definition=function(conn,...) TRUE)
#######     end kludges   ######

setClass("TSzipConnection", contains=c("zipConnection", "conType","TSdb"),
   slots=c(suffix="character") )

setMethod("TSconnect",   signature(q="zipConnection", dbname="missing"),
  definition=function(q, dbname, 
                suffix=c("Open","High","Low","Close","Volume","OI"), ...){ 
   #  user / password / host  for future consideration
   # may need to to have this function specific to dbname  cases as in TSsdmx
   dbname <- q@dbname
   
   new("TSzipConnection", dbname=dbname, 
        hasVintages=FALSE, hasPanels=FALSE,
	#read.csvArgs=list(...), 
	suffix=suffix) 
   } 
   )

setMethod("TSdates",
  signature(serIDs="character", con="TSzipConnection", vintage="ANY", panel="ANY"),
   definition= function(serIDs, con, vintage=NULL, panel=NULL, ... )  
{  # Indicate  dates for which data is available.
   # This requires retrieving series individually so they are not truncated.
   r <- av <- st <- en <- tb <- NULL
   for (i in 1:length(serIDs))
     {r <- try(TSget(serIDs[i], con), silent = TRUE)

      if(inherits(r, "try-error") ) {
        av <- c(av, FALSE)
	st <- append(st, list(NA))
	en <- append(en, list(NA))
	tb <- rbind(tb, NA)
	}
      else  {
        av <- c(av, TRUE)
        st <- append(st, list(tfstart(r)))
        en <- append(en, list(tfend(r)))
        tb <- rbind(tb,tffrequency(r))
        }
      }
  r <- serIDs
  attr(r, "TSdates") <- av
  attr(r, "start") <- st
  attr(r, "end")   <- en
  attr(r, "frequency")   <- tb
  class(r) <- "TSdates"
  r
} )

setMethod("TSget",     signature(serIDs="character", con="TSzipConnection"),
   definition=function(serIDs, con, TSrepresentation=options()$TSrepresentation,
       tf=NULL, start=tfstart(tf), end=tfend(tf),
       names=NULL, quote=con@suffix, ...){ 
    
   if(is.null(names)) names <- c(t(outer(serIDs, quote, paste, sep=".")))
   quote <- con@suffix %in% quote

   dir <- tempfile()
   dir.create(dir)
   on.exit(unlink(dir) )
   mat <- NULL
   
   for (i in 1:length(serIDs)){
      url <- paste(con@dbname, "/", serIDs[i], ".zip", sep="")
      file <- paste(dir, "/", serIDs[i], ".zip", sep="")

      zz <- try(download.file(url, file, quiet = TRUE, mode = "wb",
   		      cacheOK = TRUE),  silent=TRUE) 
      if(inherits(zz, "try-error") || (0 != zz)) 
       stop("download.file error, possibly could not find url ",  url,
            " or file ", file)

      zz <- try(unzip(file, overwrite = TRUE, exdir=dir))
      #zz <- try(system(paste("unzip", file, " -d ", dir)),  silent=TRUE)
      if(inherits(zz, "try-error")) stop("Could not unzip file ", file)

      file <- paste(dir, "/", serIDs[i], ".txt", sep="")
      zz <- try(read.csv(file),  silent=TRUE)
   		      #method=c("csv","tsv","tab"), perl="perl")
  #  # header=TRUE, sep=",", quote="\"", dec=".", fill=TRUE, comment.char=""  
  #  #  could use colClasses
     
      if(inherits(zz, "try-error")) 
   	    stop("Could read downloaded file ",  file, zz)
    
      zz <- as.matrix(zz)
      dates <- as.Date(zz[,1], format="%m/%d/%Y")
      zzz <- try(as.numeric(zz[,-1]),  silent=TRUE)
      if(inherits(zzz, "try-error")) 
   	    stop("Error converting  data to numeric.", data)
      
      d <- matrix(zzz, NROW(zz), NCOL(zz)-1 )
      d <- zoo::zoo(d[, quote], order.by=dates)
  
      mat <- tbind(mat,d)
      }
      
    mat <- tfwindow(mat, tf=tf, start=start, end=end)
    
    mat <- tframePlus::changeTSrepresentation(mat, TSrepresentation)

    seriesNames(mat) <- names
    desc <- paste(names, " from ", con@dbname)

    TSmeta(mat) <- new("TSmeta", serIDs=serIDs,  dbname=con@dbname, 
        hasVintages=con@hasVintages, hasPanels=con@hasPanels,
  	conType=class(con), DateStamp= Sys.time(), 
	TSdoc=paste(desc, "retrieved ", Sys.time()),
	TSdescription=desc,
	TSlabel=names,
	TSsource=con@dbname # could be better
	) 
    mat
    } 
    )


#setMethod("TSput",     signature(x="ANY", serIDs="character", con="TSzipConnection"),
#   definition= function(x, serIDs=seriesNames(data), con, ...)   
#    "TSput for TSzip connection not supported." )

setMethod("TSdescription",   signature(x="character", con="TSzipConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSdescription for TSzip connection not supported." )


setMethod("TSdoc",   signature(x="character", con="TSzipConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSdoc for TSzip connection not supported." )

setMethod("TSlabel",   signature(x="character", con="TSzipConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSlabel for TSzip connection not supported." )

setMethod("TSsource",   signature(x="character", con="TSzipConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSsource for TSzip connection not supported." )

