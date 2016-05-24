# these methods us a fame server connection (rather than local fame fake
#  connection to local or remote db).
   
####### some kludges to make this look like DBI. ######

setClass("fameServerDriver", contains=c("DBIDriver"), slots=c(Id = "character")) 

# require("DBI") for this
setOldClass("fameConnection",  prototype=structure(integer(1),
  class='fameConnection',      server=character(1), 
  host=character(1),     user=character(1), password=character(1)))

setClass("TSfameServerConnection",
    contains=c("DBIConnection","conType","TSdb","fameConnection"),
   slots=c(current="character"))
   # current is only used with vintages (faking info the db should have)

# these prevents error messages
setMethod("dbDisconnect", signature(conn="TSfameServerConnection"), 
   definition=function(conn,...){
      z <-  close(S3Part(conn))
      invisible(!inherits(z, "try-error"))})

#######     end kludges   ######


setMethod("TSconnect", signature(q="TSfameServerConnection", dbname="missing"),
  definition= function(q, dbname, 
     service = "", host = "", user = "", password = "", 
	      current=NA, ...){
   #It might be possible to leave the Fame db open, but getfame needs it closed.

   dbname <- q@dbname

   #ensure the db name does not end in .db, fameServer adds this.
   #dbname <- sub('$', '.db',sub('.db$', '', dbname))
   dbname <- sub('.db$', '', dbname)
   # if dbname is a vector (for vintages) it should have names so that
   #  dbname[vintage] can be used to subset it in TSget.
   nm <- names(dbname)
   if(is.null(nm) & (1 < length(dbname))) {
      nm <- as.character(seq(length(dbname)))
      warning("vintage names generated as sequence: ", nm)
      }
   names(dbname) <- nm

   if(!fameRunning()) fameStart(workingDB = FALSE)

   con <- fameConnection(service=service, host=host, 
             user=user, password=password, stopOnFail=TRUE)
   if(inherits(con, "try-error") )
       stop("Could not establish TSfameServerConnection.")
   for (i in seq(length(dbname))){
      Id <- try(fameDbOpen(dbname[i], connection=con, stopOnFail=TRUE))
      if(inherits(Id, "try-error") ) 
       stop("Could not establish TSfameServerConnection to ", dbname[i])
      fameDbClose(Id, closeConnection = FALSE) # this Id is not saved
      }
   
   new("TSfameServerConnection", drv="fameServer", dbname=dbname, 
	  hasVintages= (1 < length(dbname)), hasPanels=FALSE,
	  current=as.character(current), con) 
   } )


setMethod("TSdates",  
   signature(serIDs="character", con="TSfameServerConnection", vintage="ANY", panel="ANY"),
   definition= function(serIDs, con, vintage=getOption("TSvintage"), panel=NULL, ... )  
{  # Indicate  dates for which data is available.
   # This requires retrieving series individually so they are not truncated.
   if(con@hasVintages){
     if(is.null(vintage)) vintage <- "current"
     if(1 != length(vintage))
        stop("only one vintage must be specified for TSdates.")
     if("current" == vintage) vintage <- con@current
     if(is.na(vintage))
 	 stop("connection does not have a specified current vintage.")
     } 
   r <- av <- st <- en <- tb <- NULL
   for (i in 1:length(serIDs))
     {r <- try(TSget( serIDs[i], con=con), silent=TRUE) 
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

setMethod("TSget",     signature(serIDs="character", con="TSfameServerConnection"),
   definition=function(serIDs, con, TSrepresentation=getOption("TSrepresentation"),
       tf=NULL, start=tfstart(tf), end=tfend(tf), names=NULL, 
       TSdescription=FALSE, TSdoc=FALSE, TSlabel=FALSE, TSsource=TRUE,
       vintage=getOption("TSvintage"), ...)
{ # ... arguments unused
  if (is.null(TSrepresentation)) TSrepresentation <- "default"

  if ( 1 < sum(c(length(serIDs), length(vintage)) > 1))
   stop("Only one of serIDs or vintage can have length greater than 1.")

  if(con@hasVintages){
    if(is.null(names)) names <- 
         if ( length(vintage) > 1 ) vintage  else  serIDs 
    if(is.null(vintage)) vintage <- "current"
    vintage["current" == vintage] <- con@current
    if(any(is.na(vintage)))
        stop("connection does not have a specified current vintage.")
    # if vintage is a vector then serIDs needs to be expanded, 
    #  otherwise dbname needs to be expanded.
    dbname <- con@dbname[vintage]
    if ( 1 == length(vintage)) dbname  <- rep(dbname, length(serIDs))
    else  serIDs <- rep(serIDs, length(vintage))
    } 
  else {
    if(is.null(names)) names <- serIDs 
    dbname <- rep(con@dbname, length(serIDs))   
  }

  mat <- desc <- doc <- label <- source <-  rp <- NULL
  for (i in seq(length(serIDs))) {
    r <- getfame(serIDs[i], db=dbname[i], connection=S3Part(con),
              save = FALSE, envir = parent.frame(),
             start = NULL, end = NULL, getDoc = FALSE)
    if(0==length(r))
       stop("Fame retrieval failed. Series may not exist on ",dbname[i],".")
    # r is class tis
#    r <-  if((TSrepresentation=="default" | TSrepresentation=="ts")
#             && frequency(r) %in% c(1,4,12,2)) as.ts(r[[1]]) else as.zoo(r[[1]])

    if(TSrepresentation=="tis") r <- r[[1]]
    else if((TSrepresentation=="default" | TSrepresentation=="ts")
             && tif(r[[1]]) %in% c(1044,1027,1032,1050)) r <-  as.ts(r[[1]]) 
    else {
       rp <- c(rp, tifName(r[[1]]))
       r <- zoo(c(r[[1]]), order.by=as.Date(ti(r[[1]])), frequency=frequency(r[[1]]))
       }
    mat <- tbind(mat, r)
    if(TSdescription) desc   <- c(desc,   TSdescription(serIDs[i],con) ) 
    if(TSdoc)         doc    <- c(doc,    TSdoc(serIDs[i],con) ) 
    if(TSlabel)       label  <- c(label,  as(NA, "character")) #TSlabel(serIDs[i],con) ) 
    if(TSsource)      source <- c(source, "Fame db") #could be better
    }

  if(TSlabel) warning("TSlabel not supported in Fame.") 
  if (NCOL(mat) != length(serIDs)) stop("Error retrieving series", serIDs) 

  mat <- tfwindow(mat, tf=tf, start=start, end=end)

  if( (!is.null(rp)) && !all(is.na(rp)) ) TSrefperiod(mat) <- rp      

  if (! TSrepresentation  %in% c( "zoo", "default", "tis")){
      require("tframePlus")
      mat <- changeTSrepresentation(mat, TSrepresentation)
      }

  seriesNames(mat) <- names  

  TSmeta(mat) <- new("TSmeta", serIDs=serIDs, dbname=dbname, 
      hasVintages=con@hasVintages, hasPanels=con@hasPanels,
      conType=class(con), 
      DateStamp=Sys.time(), 
      TSdescription=if(TSdescription) paste(desc, " from ", dbname) else NA, 
      TSdoc=if(TSdoc)        doc   else NA,
      TSlabel=if(TSlabel)   label  else NA,
      TSsource=if(TSsource) source else NA )
  mat
} )


setMethod("TSput",     signature(x="ANY", serIDs="character", con="TSfameServerConnection"),
   definition= function(x, serIDs=seriesNames(x), con,   
       TSdescription.=TSdescription(x), TSdoc.=TSdoc(x), TSlabel.=NULL, 
       TSsource.=NULL, warn=TRUE, ...) 
 {
  if (con@hasVintages)
    stop("TSput does not support vintages. Open the con to a single dbname.")
  if (!is.null(TSlabel.))  warning("TSlabel is not supported in Fame.")
  if (!is.null(TSsource.)) warning("TSsource is not supported in Fame.")
  if ( 1 < length(con@dbname))
    stop("TSput does not support vintages. Open the con to a single dbname.")
  ids <-  serIDs 
  x <- as.matrix(as.tis(x)) # clobbers seriesNames(x)
  #ids <- gsub(" ", "", serIDs ) # remove spaces in id
  #if(! all( ids == serIDs)) warning("spaces removed from series names.")
  #rP <- TSrefperiod(x)
  #N <- periods(x)
  ok <- TRUE
  for (i in ids)  ok <- ok & !TSexists(i, con=con)

  if (warn & !ok) warning("error series already exist on database.")
  
  if (ok) {
    Id <- fameDbOpen(con@dbname, accessMode = "update")
    on.exit(fameDbClose(Id))
    if (ok) for (i in seq(length(ids))) {
      v <- x[,i]
      documentation(v) <- TSdoc.[i]
      description(v) <- TSdescription.[i]
      #putfame does not write doc and des. fameWriteSeries needs open and close
      ok <- ok & 0==fameWriteSeries(Id, ids[i], v,
  		       update=FALSE, checkBasisAndObserved=FALSE)
      }
  }
  
  if (warn & !ok) warning("error putting data on database.")
  new("logicalId",  ok, 
       TSid=new("TSid", serIDs=serIDs, dbname=con@dbname, 
         conType=class(con), hasVintages=con@hasVintages, hasPanels=con@hasPanels,
	 DateStamp=Sys.time()))
  } )



setMethod("TSdescription",   signature(x="character", con="TSfameServerConnection"),
   definition= function(x, con=getOption("TSconnection"), ...){
     r <- fameWhats(con@dbname[1], x, connection=S3Part(con), getDoc = TRUE)$des 
     if (is.null(r)) stop("Series (probably) does not exist.")
     #if(is.null(r) || is.na(r)|| ("NA" == r)) NA else r 
     if(is.na(r)|| ("NA" == r)) NA else r })

setMethod("TSdoc",   signature(x="character", con="TSfameServerConnection"),
   definition= function(x, con=getOption("TSconnection"), ...){
     r <- fameWhats(con@dbname[1], x, connection=S3Part(con), getDoc = TRUE)$doc
     if (is.null(r)) stop("Series (probably) does not exist.")
     #if(is.null(r) || is.na(r)|| ("NA" == r)) NA else r 
     if(is.na(r)|| ("NA" == r)) NA else r })

#TSlabel,TSsource, get used for new("Meta", so issuing a warning is not a good idea here.

setMethod("TSlabel",   signature(x="character", con="TSfameServerConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
     as(NA, "character") )

setMethod("TSsource",   signature(x="character", con="TSfameServerConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
     as(NA, "character") )

setMethod("TSdelete",
   signature(serIDs="character", con="TSfameServerConnection", vintage="ANY", panel="ANY"),
   definition= function(serIDs, con=getOption("TSconnection"),  
            vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...){
   if (con@hasVintages)
       stop("TSdelete does not support vintages. Open the con to a single dbname.")
   Id <- try(fameDbOpen(con@dbname, connection=S3Part(con), stopOnFail=TRUE))
   if(inherits(Id, "try-error") )
       stop("Could not establish TSfameServerConnection to ", con@dbname)
   
    ok <- TRUE
    for (i in seq(length(serIDs))) 
      ok <- ok & 0 == fameDeleteObject(Id, serIDs[i]) 

    fameDbClose(Id, closeConnection = FALSE) # Id is not saved

    new("logicalId",  ok, 
         TSid=new("TSid", serIDs=serIDs, dbname=con@dbname, 
           conType=class(con), hasVintages=con@hasVintages, hasPanels=con@hasPanels,
	   DateStamp=Sys.time()))
   })


setMethod("TSexists", 
 signature(serIDs="character", con="TSfameServerConnection", vintage="ANY", panel="ANY"),
 definition= function(serIDs, con=getOption("TSconnection"), 
                      vintage=NULL, panel=NULL, ...){
   if (con@hasVintages)
       stop("TSexists does not support vintages. Open the con to a single dbname.")
   op <- options(warn=-1)
   on.exit(options(op))
   ok <- fameWhats(con@dbname, serIDs, connection=S3Part(con), getDoc = FALSE)
   new("logicalId",  !is.null(ok), 
       TSid=new("TSid", serIDs=serIDs, dbname=con@dbname, 
         conType=class(con), hasVintages=con@hasVintages, hasPanels=con@hasPanels,
	 DateStamp=NA))
   })

setMethod("TSvintages",  
   signature(con="TSfameServerConnection"), 
   definition= function(con){
     if(!con@hasVintages) NULL else names(con@dbname)
     } ) 

