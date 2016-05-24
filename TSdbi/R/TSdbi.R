#setClassUnion("OptionalPOSIXct",   c("POSIXct",   "NULL"))
# bug work around
setClassUnion("OptionalPOSIXct",   c("POSIXct",   "logical"))


# inherit this into TSconnections so methods can abstract  from
#   specific drivers 
setClass("conType", slots = c( q="character"), contains="VIRTUAL" )

# this just has db info
setClass("TSdb", slots = c( dbname="character", 
    hasVintages="logical", hasPanels="logical"),  contains="VIRTUAL" )

#The connection type is the class of an actual connection, which includes the
# above virtual class. The conType in next TSid duplicates this, but otherwise
# I don't see how to put the info in the TSid.

# this has serIDs as well as there source db (so could be used to retrieve data)
setClass("TSid",  contains="TSdb", slots = c(serIDs="character", 
                     conType="character",DateStamp="OptionalPOSIXct")) 

# TSmeta has serIDs and source db (in TSid) as well as documentation
setClassUnion("OptionalChar",   c("character",   "logical"))
#  NA is logical so this allows NA
setClass("TSmeta", contains="TSid", 
  slots = c(TSdescription="OptionalChar",    TSdoc="OptionalChar",
                  TSlabel="OptionalChar", TSsource="OptionalChar")) 

# this is a logical but has TSid (so could be used to retrieve data)
setClass("logicalId", contains="logical", slots = c(TSid="TSid")) 

setGeneric("TSmeta", def= function(x, con=getOption("TSconnection"), ...)
    standardGeneric("TSmeta"))

setMethod("show", "logicalId", function(object) show(object@.Data))


# this really should be for con="TSconnection" rather than ANY
setMethod("TSmeta",   signature(x="character", con="ANY"),
   definition= function(x, con=getOption("TSconnection"), ...){
       new("TSmeta", serIDs=x, dbname=con@dbname, 
               hasVintages=con@hasVintages, hasPanels=con@hasPanels,
               conType=class(con), DateStamp=NA, 
	       TSdoc=TSdoc(x, con=con, ...),
               TSdescription=TSdescription(x, con=con, ...),
	       TSlabel=TSlabel(x, con=con, ...),
	       TSsource=TSsource(x, con=con, ...))
	} )

setMethod("TSmeta",   signature(x="character", con="missing"),
   definition= function(x, con=getOption("TSconnection"), ...) 
       TSmeta(x, con=getOption("TSconnection"), ...) )

# extract meta data from an object, not from the db
setMethod("TSmeta",   signature(x="ANY", con="missing"),
   definition=  function(x, con, ...) {
      m <- attr(x, "TSmeta")
      if(is.null(m)) new("TSmeta",serIDs=seriesNames(x), 
	         dbname="", hasVintages=FALSE, hasPanels=FALSE,
	         conType="", 
		 DateStamp     = NA, 
		 TSdescription = NA, 
		 TSdoc         = NA,  
		 TSlabel       = NA,
		 TSsource      = NA) else m
      })

setGeneric("TSmeta<-", 
   def= function(x, value) standardGeneric("TSmeta<-"),
   useAsDefault= function(x, value){
      if (!is(value,"TSmeta")) 
         stop("trying to set attribute TSmeta incorrectly.") 
      attr(x, "TSmeta") <- value
      x
      })

#setMethod("show", "TSdb", function(object) {
#    cat("database ", object@dbname) 
#    if (object@vintage) cat( " Has vintages." )
#    if (object@panel) cat( " Has panels." )
#    cat("\n") 
#    invisible(object)
#    })
#
#setMethod("print", "TSdb", function(x, ...) {
#    cat("database ", x@dbname) 
#    if (x@vintage) cat( " Has vintages." )
#    if (x@panel) cat( " Has panels." )
#    cat("\n") 
#    invisible(x)
#    })

setGeneric("TSconnect", def= function(q, dbname, ...) standardGeneric("TSconnect"))

# note that ... is passed to the driver, dbConnect, and TSconnect in next
setMethod("TSconnect",   signature(q="character", dbname="character"),
   definition=function(q, dbname, ...) TSconnect(dbConnect(
        get(q, mode="function")(...), dbname=dbname, ... ), ...))

# This would work above for DBI sql engines, but for others the TS* package needs to export
#        TSconnect(dbConnect(getExportedValue(paste("R", q,sep=""), q)(),
#	           dbname=dbname) ))


######### TSdescription #########

setGeneric("TSdescription<-", 
   def= function(x, value) standardGeneric("TSdescription<-"),
   useAsDefault=  function (x, value){
    m <- TSmeta(x)
    m@TSdescription <- value
    TSmeta(x) <- m
    x})

setGeneric("TSdescription", def= function(x, con=getOption("TSconnection"), ...)
    standardGeneric("TSdescription"))

setMethod("TSdescription",   signature(x="character", con="missing"),
   definition= function(x, con=getOption("TSconnection"), ...) 
       TSdescription(x, con=getOption("TSconnection"), ...) )

setMethod("TSdescription",   signature(x="ANY", con="missing"),
   definition=  function(x, con, ...)TSmeta(x)@TSdescription)

#  next is for case where there is no method for con  
setMethod("TSdescription",   signature(x="character", con="ANY"),
   definition= function(x, con=getOption("TSconnection"), ...) {
       if(is.null(con)) stop("NULL con is not allowed. See ?TSdescription.")
       else stop("con class ", class(con), 
       " is not supported. (Check this is a TSdbi connection, not a raw DBI connection.)")} )

#  Next two methods are for case where the user mistakenly specifies
#     serIDS="whatever"  rather than x="whatever"
#   ( A natural mistaken as this is the syntax for other functions.)
setMethod("TSdescription",   signature(x="missing", con="ANY"),
   definition=  function(x, con, serIDs, ...) {
        if(missing(serIDs)) stop("missing argument x, must be specified.")
	else TSdescription(x=serIDs, con=con, ...)
	})

setMethod("TSdescription",   signature(x="missing", con="missing"),
   definition=  function(x, serIDs, ...) {
        if(missing(serIDs)) stop("missing argument x, must be specified.")
	else TSdescription(x=serIDs, ...)
	})

######### TSdoc #########

setGeneric("TSdoc<-", 
   def= function(x, value) standardGeneric("TSdoc<-"),
   useAsDefault=  function (x, value){
    m <- TSmeta(x)
    m@TSdoc <- value
    TSmeta(x) <- m
    x})

setGeneric("TSdoc", 
   def= function(x, con=getOption("TSconnection"), ...) standardGeneric("TSdoc"))

setMethod("TSdoc",   signature(x="character", con="missing"),
   definition= function(x, con=getOption("TSconnection"), ...){ 
       if(is.null(con)) 
	  stop("con should be specified or set with options(TSconnection=con). See ?TSdoc.") 
       TSdoc(x, con=con, ...)} )

setMethod("TSdoc",   signature(x="ANY", con="missing"),
   definition=  function(x, con, ...) TSmeta(x)@TSdoc)

#  next is for case where there is no method for con  
setMethod("TSdoc",   signature(x="character", con="ANY"),
   definition= function(x, con=getOption("TSconnection"), ...) {
       if(is.null(con)) stop("NULL con is not allowed. See ?TSdoc.")
       else stop("con class ", class(con), 
       " is not supported. (Check this is a TSdbi connection, not a raw DBI connection.)")} )

#  Next two methods are for case where the user mistakenly specifies
#     serIDS="whatever"  rather than x="whatever"
#   ( A natural mistaken as this is the syntax for other functions.)
setMethod("TSdoc",   signature(x="missing", con="ANY"),
   definition=  function(x, con, serIDs, ...) {
        if(missing(serIDs)) stop("missing argument x, must be specified.")
	else TSdoc(x=serIDs, con=con, ...)
	})

setMethod("TSdoc",   signature(x="missing", con="missing"),
   definition=  function(x, serIDs, ...) {
        if(missing(serIDs)) stop("missing argument x, must be specified.")
	else TSdoc(x=serIDs, ...)
	})

######### TSlabel #########

setGeneric("TSlabel<-", 
   def= function(x, value) standardGeneric("TSlabel<-"),
   useAsDefault=  function (x, value){
    m <- TSmeta(x)
    m@TSlabel <- value
    TSmeta(x) <- m
    x})

setGeneric("TSlabel", 
   def= function(x, con=getOption("TSconnection"), ...) standardGeneric("TSlabel"))

setMethod("TSlabel",   signature(x="character", con="missing"),
   definition= function(x, con=getOption("TSconnection"), ...){ 
       if(is.null(con)) 
          stop("con should be specified or set with options(TSconnection=con). See ?TSlabel.") 
       TSlabel(x, con=con, ...)} )

setMethod("TSlabel",   signature(x="ANY", con="missing"),
   definition=  function(x, con, ...) TSmeta(x)@TSlabel)

#  next is for case where there is no method for con  
setMethod("TSlabel",   signature(x="character", con="ANY"),
   definition= function(x, con=getOption("TSconnection"), ...) {
       if(is.null(con)) stop("NULL con is not allowed. See ?TSlabel.")
       else stop("con class ", class(con), 
       " is not supported. (Check this is a TSdbi connection, not a raw DBI connection.)")} )

#  Next two methods are for case where the user mistakenly specifies
#     serIDS="whatever"  rather than x="whatever"
#   ( A natural mistaken as this is the syntax for other functions.)
setMethod("TSlabel",   signature(x="missing", con="ANY"),
   definition=  function(x, con, serIDs, ...) {
        if(missing(serIDs)) stop("missing argument x, must be specified.")
	else TSlabel(x=serIDs, con=con, ...)
	})

setMethod("TSlabel",   signature(x="missing", con="missing"),
   definition=  function(x, serIDs, ...) {
        if(missing(serIDs)) stop("missing argument x, must be specified.")
	else TSlabel(x=serIDs, ...)
	})

######### TSsource #########

setGeneric("TSsource<-", 
   def= function(x, value) standardGeneric("TSsource<-"),
   useAsDefault=  function (x, value){
    m <- TSmeta(x)
    m@TSsource <- value
    TSmeta(x) <- m
    x})


setGeneric("TSsource", 
   def= function(x, con=getOption("TSconnection"), ...) standardGeneric("TSsource"))

setMethod("TSsource",   signature(x="character", con="missing"),
   definition= function(x, con=getOption("TSconnection"), ...){ 
       if(is.null(con)) 
          stop("con should be specified or set with options(TSconnection=con). See ?TSsource.") 
       TSsource(x, con=con, ...)} )

setMethod("TSsource",   signature(x="ANY", con="missing"),
   definition=  function(x, con, ...) TSmeta(x)@TSsource)

#  next is for case where there is no method for con  
setMethod("TSsource",   signature(x="character", con="ANY"),
   definition= function(x, con=getOption("TSconnection"), ...) {
       if(is.null(con)) stop("NULL con is not allowed. See ?TSsource.")
       else stop("con class ", class(con), 
       " is not supported. (Check this is a TSdbi connection, not a raw DBI connection.)")} )

#  Next two methods are for case where the user mistakenly specifies
#     serIDS="whatever"  rather than x="whatever"
#   ( A natural mistaken as this is the syntax for other functions.)
setMethod("TSsource",   signature(x="missing", con="ANY"),
   definition=  function(x, con, serIDs, ...) {
        if(missing(serIDs)) stop("missing argument x, must be specified.")
	else TSsource(x=serIDs, con=con, ...)
	})

setMethod("TSsource",   signature(x="missing", con="missing"),
   definition=  function(x, serIDs, ...) {
        if(missing(serIDs)) stop("missing argument x, must be specified.")
	else TSsource(x=serIDs, ...)
	})

setMethod("show", "TSmeta", function(object) {
    cat("serIDs: ", object@serIDs,"\n") 
    if(!all(is.na(object@TSsource))) cat("source: ", object@TSsource,"\n") 
    cat(" from dbname ", object@dbname) 
    cat(" using ", object@conType) 
        if(!is.na(object@DateStamp)) cat(" on ",as.character(object@DateStamp)) 
    cat("\n") 
    if(!all(is.na(object@TSlabel)))       
            cat("label:        ", object@TSlabel, "\n") 
    if(!all(is.na(object@TSdescription)))
            cat("description: ", object@TSdescription, "\n") 
    if(!all(is.na(object@TSdoc)) )
            cat("documentaion: ", object@TSdoc, "\n") 
    invisible(NULL)
    })

#setMethod("print", "TSmeta", function(x, ...) {
#    cat("serIDs: ", x@serIDs, " from dbname ", x@dbname) 
#    cat("(type: ", x@conType, ")\n") 
#    cat("description: ", x@TSdescription, "\n") 
#    cat("documentaion: ", x@TSdoc, "\n") 
#    invisible(x)
#    })


setGeneric("TSrefperiod", 
   def= function(x) standardGeneric("TSrefperiod"),
   useAsDefault= function(x) attr(x, "TSrefperiod"))

setGeneric("TSrefperiod<-", 
   def= function(x, value) standardGeneric("TSrefperiod<-"),
   useAsDefault= function (x, value){attr(x, "TSrefperiod") <- value ; x })


TSseriesIDs      <- function(x)TSmeta(x)@serIDs
TScon            <- function(x)TSmeta(x)@con
TSextractionDate <- function(x)TSmeta(x)@DateStamp

# host = "" is interpreted as localhost by the driver (see ?mysqlNewConnection)
#see  vignette("zoo") for other examples


setGeneric("TSput", valueClass="logicalId",   
   def= function(x, serIDs=seriesNames(x), con=getOption("TSconnection"), ...)
             standardGeneric("TSput")) 

setMethod("TSput",   signature(x="ANY", serIDs="character", con="missing"),
   definition= function(x, serIDs, con=getOption("TSconnection"), ...){ 
       if(is.null(con)) 
          stop("con should be specified or set with options(TSconnection=con). See ?TSput.") 
       TSput(x, serIDs=serIDs, con=con, ...)})

setMethod("TSput",   signature(x="ANY", serIDs="missing", con="missing"),
   definition= function(x, serIDs, con=getOption("TSconnection"), ...){ 
       if(is.null(con)) 
          stop("con should be specified or set with options(TSconnection=con). See ?TSput.")  
       TSput(x, serIDs=seriesNames(x), con=con, ...)})

setMethod("TSput",   signature(x="ANY", serIDs="DBIConnection", con="missing"),
   definition= function(x, serIDs, con, ...)
          TSput(x, serIDs=seriesNames(x), con=serIDs, ...))

#  next is for case where there is no method for con  
setMethod("TSput",   signature(x="ANY", serIDs="character", con="ANY"),
   definition= function(x, serIDs=seriesNames(x), con=getOption("TSconnection"), ...) {
       if(is.null(con)) stop("NULL con is not allowed. See ?TSput.")
       else stop("con class ", class(con), 
       " is not supported. (Check this is a TSdbi connection, not a raw DBI connection.)")} )

   
setGeneric("TSreplace",  valueClass="logicalId",
   def= function(x, serIDs=seriesNames(x), con=getOption("TSconnection"), 
           vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
              standardGeneric("TSreplace"),
   useAsDefault= function(x,serIDs=seriesNames(x), con=getOption("TSconnection"), 
           vintage=getOption("TSvintage"), panel=getOption("TSpanel"),  ...) {
      if(missing(con) & (!missing(serIDs)) && is(serIDs, "DBIConnection")) 
	     return(TSreplace(x, serIDs=seriesNames(x), con=serIDs,
	             vintage=vintage, panel=panel, ...))
      for (s in serIDs) {
        if(TSexists(serIDs=s, con=con, vintage=vintage, panel=panel,...))
           TSdelete(serIDs=s, con=con, vintage=vintage, panel=panel, ...)
	} 
      TSput(x, serIDs=serIDs, con=con, vintage=vintage, panel=panel, ...)
   })


setGeneric("TSdelete",  valueClass="logicalId",
   def= function(serIDs, con=getOption("TSconnection"), 
           vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
              standardGeneric("TSdelete") )

setMethod("TSdelete",   
   signature(serIDs="character", con="missing", vintage="ANY", panel="ANY"),
   definition= function(serIDs, con=getOption("TSconnection"), 
           vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...){ 
       if(is.null(con)) 
          stop("con should be specified or set with options(TSconnection=con). See ?TSdelete.") 
       TSdelete(serIDs, con=con, vintage=vintage, panel=panel, ...)} )

#  next is for case where there is no method for con  
setMethod("TSdelete",   
   signature(serIDs="character", con="ANY", vintage="ANY", panel="ANY"),
   definition= function(serIDs, con=getOption("TSconnection"),
        vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...) {
       if(is.null(con)) stop("NULL con is not allowed. See ?TSdelete.")
       else stop("con class ", class(con), 
       " is not supported. (Check this is a TSdbi connection, not a raw DBI connection.)")} )

   

setGeneric("TSget", 
   def= function(serIDs, con=getOption("TSconnection"), ...)
           standardGeneric("TSget") )

setMethod("TSget",   signature(serIDs="character", con="missing"),
   definition= function(serIDs, con=getOption("TSconnection"), ...){ 
       if(is.null(con)) 
          stop("con should be specified or set with options(TSconnection=con). See ?TSget.")
       TSget(serIDs, con=con, ...)} )

#  next is for case where there is no method for con  
setMethod("TSget",   signature(serIDs="character", con="ANY"),
   definition= function(serIDs, con=getOption("TSconnection"), ...){
       if(is.null(con)) stop("NULL con is not allowed. See ?TSget.")
       else stop("con class ", class(con), 
       " is not supported. (Check this is a TSdbi connection, not a raw DBI connection.)")} )

setGeneric("TSdates", def=
   function(serIDs, con=getOption("TSconnection"), 
           vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
   standardGeneric("TSdates") )

setMethod("TSdates",
   signature(serIDs="character", con="missing", vintage="ANY", panel="ANY"),
   definition= function(serIDs, con=getOption("TSconnection"), 
           vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...){ 
       if(is.null(con)) 
          stop("con should be specified or set with options(TSconnection=con). See ?TSdates.") 
       TSdates(serIDs, con=con, vintage=vintage, panel=panel, ...)} )

#  next is for case where there is no method for con  
setMethod("TSdates",
   signature(serIDs="character", con="ANY", vintage="ANY", panel="ANY"),
   definition= function(serIDs, con=getOption("TSconnection"), 
           vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...) {
       if(is.null(con)) stop("NULL con is not allowed. See ?TSdates.")
       else stop("con class ", class(con), 
       " is not supported. (Check this is a TSdbi connection, not a raw DBI connection.)")} )

start.TSdates   <- function(x, ...) attr(x, "start")
tfstart.TSdates <- function(x) attr(x, "start")
end.TSdates   <- function(x, ...) attr(x, "end")
tfend.TSdates <- function(x) attr(x, "end")

print.TSdates   <- function(x, ...) {
   r <- NULL
   av <-  attr(x, "TSdates")    # a logical indicator of availability
   st <-  attr(x, "start")
   en <-  attr(x, "end")  
   fr <-  attr(x, "frequency")  # used by defunct TSpadi, possibly elsewhere
   tb <-  attr(x, "tbl")  
   rP <-  attr(x, "TSrefperiod")
   for (i in 1:length(x)) {
      if (!av[i]) r <- c(r, paste(x[i], " not available"))
      else {
        extra <- paste(tb[i], " ", rP[i], " ", fr[i])
        r <- c(r, paste(x[i], "from", paste(st[[i]], collapse=" "),
	                        "to", paste(en[[i]], collapse=" "), extra))
	}
      }
   print(as.matrix(r), ...)
   }
     

setGeneric("TSexists", valueClass="logicalId",
 def= function(serIDs, con=getOption("TSconnection"), 
           vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
 	    standardGeneric("TSexists"),
 useAsDefault=function(serIDs, con=getOption("TSconnection"), 
           vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...){
    #vintages treated as a synonym for vintage.
    #This is done below using TSdates, so the handling of SQL (or not) is looked
    # after in TSdates. It could be made into a direct query like:
    #    vintage %in% .. SELECT vintage FROM Meta where id = "V122707";
    # but that would require a method in each package as for TSdates.
    
    # catch a possible common mistake
    if("vintages" %in% names(list(...)))
       stop("TSexists argument is 'vintage' not 'vintages'.")
    
    if(is.null(vintage)) {
      z <-  try(TSdates(serIDs,
               con, vintage=NULL, panel=panel, ...), silent=TRUE)
      av <- if( inherits(z, "try-error")) FALSE else all(attr(z, "TSdates"))
      }
    else {
      av <- NULL
      for (i in 1:length(vintage)) {
         z <- try(TSdates(serIDs, 
	       con, vintage=vintage[i], panel=panel, ...), silent=TRUE)
         av <- c(av, 
	    if( inherits(z, "try-error")) FALSE else all(attr(z, "TSdates")))
         }
      }

    new("logicalId", av, 
       TSid=new("TSid", serIDs=serIDs, dbname=con@dbname, 
              conType=class(con), hasVintages=con@hasVintages, hasPanels=con@hasPanels,
	      DateStamp=NA))})


########## little utilities #######
# PostgreSQL and MySQL" failing with cannot coerce type 'S4' to vector of type 'integer'

TSfinddb <- function(dbname=NULL, driverOrder=c("MySQL", "SQLite", "PostgreSQL")) {
  if(is.null(dbname)) stop("dbname must be supplied.")
  org <- search()
  ok <- FALSE
  for (s in driverOrder) {
     pkg <- paste("TS", s,sep="")
     if(require(pkg, quietly = TRUE, character.only = TRUE)) {
        con <- try(TSconnect(s, dbname), silent = TRUE)
        ok <- ! inherits(con, "try-error")
	}
     if (ok) break
     else {
        z <- paste("package:TS", s, sep="")
	if(! z %in% org) detach(z, character.only = TRUE)
        z <- paste("package:R", s, sep="")
	if(! z %in% org) detach(z, character.only = TRUE)
        }
     } 
  if( inherits(con, "try-error")) {
      warning(dbname, " not found.")
      return(NULL)
      }
  else return(con)
  }


# return all vintage ids

setGeneric("TSvintages", 
   def=function(con=getOption("TSconnection")) standardGeneric("TSvintages"))

# setting the default this way would mean the method does not need to be
# defined in SQL packages, but requires other packages to define a method,
# otherwise they fail with a non-intuitive message.
#setGeneric("TSvintages", 
#   def=function(con=getOption("TSconnection")) standardGeneric("TSvintages"),
#   useAsDefault= function(con) {
#   if(is.null(con)) stop("NULL con is not allowed. See ?vintages.")
#   if(!con@hasVintages) return(NULL)   
#   dbGetQuery(con,"SELECT  DISTINCT(vintage) FROM  vintages;" )$vintage
#  } )

setMethod("TSvintages",
   signature(con="missing"),
   definition= function(con=getOption("TSconnection")){ 
       if(is.null(con)) 
          stop("con should be specified or set with options(TSconnection=con). See ?TSvintages.") 
       TSvintages(con=con)} )

#  next is for case where there is no method for con  
setMethod("TSvintages",
   signature(con="ANY"),
   definition= function(con=getOption("TSconnection")) {
       if(is.null(con)) stop("NULL con is not allowed. See ?TSvintages.")
       else stop("con class ", class(con), 
       " is not supported. (Check this is a TSdbi connection, not a raw DBI connection.)")} )


