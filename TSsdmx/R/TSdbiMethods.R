sdmx <- function(...) {
  drv <- "sdmx"
  attr(drv, "package") <- "TSsdmx"
  new("sdmxDriver")
  }

# there is an SDMX primer at
# http://www.ecb.int/stats/services/sdmx/html/index.en.html

# NB firebug shows browser requests to server, so is useful for seeing what is
#  sent to the server

####### some kludges to make this look like DBI. ######
#for this require("DBI") ; require("RJSDMX")

setClass("sdmxDriver", contains=c("DBIDriver")) 

setClass("sdmxConnection", contains=c("DBIConnection", "sdmxDriver"),
   slots=c(dbname="character") )

# slot getTimeSeries="function" could be here but J() function 
#     value may not persist in method across sessions??
setMethod("dbConnect", signature(drv="sdmxDriver"), 
     definition=function(drv, dbname, ...) 
         new("sdmxConnection", dbname=dbname))

# this does nothing but prevent errors if it is called. 
setMethod("dbDisconnect", signature(conn="sdmxConnection"), 
     definition=function(conn,...) TRUE)

#######     end kludges   ######
# require(TSdbi)

setClass("TSsdmxConnection", contains=c("DBIConnection", "conType","TSdb"),
   slots=c( getTimeSeries="function",
      user="character", password="character", host="character") )

setMethod("TSconnect",   signature(q="sdmxConnection", dbname="missing"),
  definition= function(q, dbname, user="", password="", host="", ...){
   #  user / password / host  for future consideration
   
   dbname <- q@dbname #  eg "BIS" "ILO" "ECB" "OECD" "EUROSTAT"

   providers <- try(getProviders())
   
   if(inherits(providers, "try-error")) 
         stop("Error trying to verify provider ",  dbname)

   if( ! (dbname %in% providers))
          stop("dbname ", dbname, " not a recognized SDMX provider.")
   
   new("TSsdmxConnection", dbname=dbname, getTimeSeries=get,
        hasVintages=FALSE, hasPanels=FALSE, 
	user=user, password=password, host=host ) 
   } )

setMethod("TSdates",
  signature(serIDs="character", con="TSsdmxConnection", vintage="ANY", panel="ANY"),
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

setMethod("TSget",     signature(serIDs="character", con="TSsdmxConnection"),
   definition  = function(serIDs, con, TSrepresentation=options()$TSrepresentation,
       tf=NULL, start=tfstart(tf), end=tfend(tf),
       names=serIDs, quiet=TRUE, ...){ 

    if (is.null(TSrepresentation)) TSrepresentation <- "ts"
    dbname <- con@dbname

    mat  <- NULL
    nm   <- NULL
    desc <- NULL
    doc  <- NULL

    if (is.character(start) & is.character(end)) {
       stenSDMX <- TRUE  # assume char is as per needed by RJSDMX
       st <- start 
       en <- end
       }
    else { 
       stenSDMX <- FALSE  # assume R dates, retrieve all and truncate below
       st <- "" 
       en <- "" 
       }
    
    #serIDs can be a vector, and a single element may have a wildcard 
    #  which returns multiple series.

    for (j in 1:length(serIDs)) {
    	ser <- try(getSDMX(dbname, serIDs[j], start=st, end=en), silent=TRUE)
       
    	if(inherits(ser, "try-error")){
    	   if(grepl("does not exist in provider", attr(ser,"condition")$message))
    	     stop(serIDs[j], " does not exist on ", dbname)
    	   else
    	     stop(serIDs[j], " error: ", attr(ser,"condition")$message)
    	   }

    	if (0 == length(ser)) stop("unknown error getting ", serIDs[j], " from ", dbname,
    		" try getSDMX('",dbname, "', '", serIDs[j], "')" )
     
    	for (i in 1:length(ser)) {
    	   mat  <- tbind(mat, ser[[i]])
    	   desc <- c(desc, names(ser[[i]]))
    	   doc  <- c(doc, attr(ser[[i]], "TITLE_COMPL"))
    	   #sdmxMeta <- append(sdmxMeta, list(attributes(ser[[i]])))
    	  }
    	desc <- c(desc, names(ser)) # one of these could be different
        nm   <- c(nm,   names(ser))
        }

    if (all(is.nan(mat))) warning("Data is all NaN.")
   
  
    mat <- tframePlus::changeTSrepresentation(mat, TSrepresentation)
    # BUG above and below usually in reverse order but this works around
    #  window BUG. May be fixed by new zoo not yet on CRAN.
    if(!stenSDMX)  mat <- tfwindow(mat, tf=tf, start=start, end=end)    
    
    # If serIDs has a * replacement skip because order is not determined.
    # If length of names does not match number of series then skip.
    # If serIDs is longer than 1 and wildcard used, skip name replacement.
    # Otherwise, for +| searches, if specified names is the right length,
    # then sort to ensure proper order.     

    wild <- any(grepl('[\\|+\\*]',serIDs))

    if (grepl('\\*',serIDs) || (length(names) != nseries(mat))) names <- nm
    else if (wild) {
       if (1 < length(serIDs)) names <- nm
       else  { # must be '[|+]'
	   parsed <- strsplit(serIDs, "\\.")[[1]]
	   m <- (1:length(parsed))[grepl("[+|]", parsed)]
	   ids <- strsplit(parsed[m], '[+|]')[[1]]
	   pre  <- if (m > 1) paste(parsed[1:(m-1)], collapse=".") else NULL
	   post <- if (m < length(parsed)) 
	            paste(parsed[(m+1):length(parsed)], collapse=".") else NULL
	   ids <- paste(pre, ids, post, sep=".")
	   if (!all(ids %in% nm)){
	       warning("skipping renaming. something is not working.")
	       names <- nm
	       }
           else {
	       dimnames(mat) <- list(NULL,nm)
	       mat <- mat[, ids] # reoder
	       }
	   }
       }

    seriesNames(mat) <- names

    TSmeta(mat) <- new("TSmeta", serIDs=serIDs,  dbname=con@dbname, 
        hasVintages=con@hasVintages, hasPanels=con@hasPanels,
  	conType=class(con), DateStamp= Sys.time(), 
	#TSdoc=paste(desc, " from ", con@dbname, "retrieved ", Sys.time()),
	TSdoc=if (is.null(doc)) "" else doc , #sdmxMeta,
	TSdescription=paste(desc, " from ", dbname),
	TSlabel=desc, 
	TSsource=con@dbname # could be better
	) 
    mat
    } 
    )


# setMethod("TSput",  signature(x="ANY", serIDs="character", con="TSsdmxConnection"),
#   definition= function(x, serIDs=seriesNames(data), con, ...)   
#    "TSput for TSsdmx connection not supported." )

setMethod("TSdescription",   signature(x="character", con="TSsdmxConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSdescription for TSsdmx connection not supported." )


setMethod("TSdoc",   signature(x="character", con="TSsdmxConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSdoc for TSsdmx connection not supported." )

setMethod("TSlabel",   signature(x="character", con="TSsdmxConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSlabel for TSsdmx connection not supported." )

setMethod("TSsource",   signature(x="character", con="TSsdmxConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        "TSsource for TSsdmx connection not supported." )


############################################################################
##############               utilities                        ##############   
############################################################################

verifyQuery <- function(provider, Q, verbose=TRUE){
  parsed <- strsplit(Q, "\\.")[[1]]
  flow <- parsed[1]
  dq   <- parsed[-1]

  if(!provider %in% getProviders()) {
     if (verbose) cat("provider", provider," not supported.\n")
     return(invisible(FALSE))
     }

  if(!flow %in% names(getFlows(provider))) {
     if (verbose) cat("flow", flow," not in providers flows.\n")
     return(invisible(FALSE))
     }

  dm <- names(getDimensions(provider, flow)) 

  if(length(dm) != length(dq)) {
     if (verbose){
         cat("query dimension not equal dimension of provider flow.\n")
         cat("query : ", dq, ".\n")
         cat("fields: ", dm, ".\n")
	 }
     return(invisible(FALSE))
     }

  for (i in 1:length(dm)) {
     z <- names(getCodes(provider, flow, dm[i]))
     if ((! (dq[i] %in% c('*','+', ''))) && (! (dq[i] %in% z))) {
           if (verbose) {
              cat("query field '", dq[i], "' not in dimension options.\n")
              cat("provider options for field ", dm[i], ": ",z,".\n")
	      }
           return(invisible(FALSE))
           }   
     }
  invisible(TRUE)
  }

#  verifyQuery('IMFx', 'PGI.CA.*.*.*.*')
#  verifyQuery('IMF', 'PGI.CA.*.*.*.*')
#  verifyQuery('IMF', 'PGI.CAN.*.*.*.*')
#  verifyQuery('IMF', 'PGI.CA..*.*.*')
  

hasData <- function(x, quiet=FALSE){
  nm <- seriesNames(x)
  N <- length(nm)
  ok <- rep(NA, N)
  for (i in 1:N){
    if(all(is.nan(x[,i]) | is.na(x[,i]))) {
       ok[i] <- FALSE
       if(!quiet) warning(nm[i], " has no data.")
       }
    else ok[i] <- TRUE
    }
   ok
   }

hasDataCount <- function(x) {
   ok <- hasData(x, quiet=TRUE)
   cat(sum(ok), " of ", length(ok), " have data.\n")
   invisible(ok)
   }

hasDataNames <- function(x) seriesNames(x)[hasData(x, quiet=TRUE)]

hasDataDescriptions <- function(x) {
  nm <- seriesNames(x)
  N <- length(nm)
  ok <- hasData(x, quiet=TRUE)

  #  this should get SDMX description
  desc <- TSmeta(x)@TSdescription
  # but this is just the flow description
  # nm <- getFlows('EUROSTAT')
  # #length(names(nm))  #5720
  # nm["ei_nama_q" == names(nm)]
  
  r <- NULL
  #for (i in 1:N) if(ok[i]) r <- rbind(r, paste(nm[i], ": ", desc[i], sep=""))
  for (i in 1:N) if(ok[i]) r <- rbind(r, desc[i])
  r
  }

