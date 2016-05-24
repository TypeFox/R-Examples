odbc <- function (...) {
   # this is really fake
   new("odbcDriver")    #DBI Id inot used
   } 

#RODBC  seems to do both driver and connection in DBI terminology
#RODBC is non-neg numeric value (or -1 for failure) with S3 class and attributes

####### some kludges to make this look like DBI. ######

#this works but the connection gives warnings when it eventually gets closed
#setOldClass("RODBC",  prototype=odbcConnect("test"))

#handle_ptr= <externalptr> 
setOldClass("RODBC",  prototype=structure(integer(1),
  class='RODBC',           connection.string=character(1), 
  handle_ptr=integer(1),   case=character(1),     id=integer(1), 
  believeNRows=logical(1), bulk_add=character(1), colQuote=character(1), 
  tabQuote=character(1),   encoding=character(1), rows_at_time=1000))

setMethod("dbListTables", signature(conn="RODBC"), 
    definition=function(conn,...) 
        as(sqlTables(channel=conn)$TABLE_NAME, "character"))

setMethod("dbExistsTable", signature(conn="RODBC", name="character"),
   definition=function(conn, name, ...) name %in% dbListTables(conn))

setMethod("dbRemoveTable", signature(conn="RODBC", name="character"),
   definition=function(conn, name, ...) {
     if (-1 == sqlDrop(conn, name, errors = FALSE) ) FALSE else TRUE})

setMethod("dbGetQuery", signature(conn="RODBC", statement="character"),
   definition = function (conn, statement, ...){
      r <- sqlQuery(channel=conn, statement, ...)
      #if( NROW(r) == 0) NULL else r
      # Above was the last line. Next is work around for odbc 
      #  returning factors rather than numeric
      if( NROW(r) == 0) return(NULL)
      if(is.null(dim(r))) return(r)  # sql other than select, or error message
      if ('tbl' %in%  names(r)) return(r) # meta query
      for (i in 1:NCOL(r)) 
       if (inherits(r[,i], "factor")) r[,i] <- as.numeric(as.character(r[,i]))
      r
      }) 

setClass("odbcDriver", contains=c("RODBC"))

setClass("odbcConnection", contains=c("RODBC", "DBIConnection"),
   slots=c(dbname="character") )

# this does not really connect, that is done in TSconnect
setMethod("dbConnect", signature(drv="odbcDriver"), 
     definition=function(drv, dbname, ...) new("odbcConnection", drv, dbname=dbname))

setMethod("dbDisconnect", signature(conn="odbcConnection"), 
    definition=function(conn,...) odbcClose(channel=conn))

#  this is pretty bad
setMethod("dbGetException", signature(conn="odbcConnection"),
   definition=function(conn, ...) list(errorNum=0))

#######     end kludges   ######

setClass("TSodbcConnection", contains=c("odbcConnection","conType","TSdb"))

setMethod("TSconnect",   signature(q="odbcConnection", dbname="missing"),
   definition=function(q, dbname, ...) {
        # q should be a connection but I cannot seem to get a valid RODBC (S3) connection
	#  passed through, so this establishes the connection here, 
	#   in contrast to other TS* packages
	nm <- q@dbname
	#nm <- as.character(attr(q "connection.string"))
        #nm <- grep("DATABASE=", scan(text=nm, sep=";", what="char", quiet=TRUE), value=TRUE)
        #nm <- gsub("DATABASE=","", nm)
 	con <- odbcConnect(dsn=nm, ...) #, uid = "", pwd = "", ...)
	if(con == -1) stop("error establishing ODBC connection.") 
	if(0 == length(dbListTables(con))){
	  odbcClose(con)
          stop("Database ",dbname," has no tables.")
	  }
	if(!dbExistsTable(con, "meta")){
	  odbcClose(con)
          stop("Database ",dbname," does not appear to be a TS database.")
	  }
	new("TSodbcConnection" , con, dbname=nm, 
  	       hasVintages=dbExistsTable(con, "vintageAlias"), 
  	       hasPanels  =dbExistsTable(con, "panels")) 
	})

setMethod("TSput",   signature(x="ANY", serIDs="character", con="TSodbcConnection"),
   definition= function(x, serIDs, con=getOption("TSconnection"), Table=NULL, 
       TSdescription.=TSdescription(x), TSdoc.=TSdoc(x), TSlabel.=TSlabel(x),
         TSsource.=TSsource(x),  
       vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
  TSputSQL(x, serIDs, con, Table=Table, 
   TSdescription.=TSdescription., TSdoc.=TSdoc., TSlabel.=TSlabel.,
     TSsource.=TSsource., 
   vintage=vintage, panel=panel) )

setMethod("TSget",   signature(serIDs="character", con="TSodbcConnection"),
   definition= function(serIDs, con=getOption("TSconnection"), 
       TSrepresentation=getOption("TSrepresentation"),
       tf=NULL, start=tfstart(tf), end=tfend(tf), names=NULL, 
       TSdescription=FALSE, TSdoc=FALSE, TSlabel=FALSE, TSsource=TRUE,
       vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
   TSgetSQL(serIDs, con, TSrepresentation=TSrepresentation,
       tf=tf, start=start, end=end,
       names=names, TSdescription=TSdescription, TSdoc=TSdoc, TSlabel=TSlabel,
         TSsource=TSsource,
       vintage=vintage, panel=panel) )

setMethod("TSdates",    signature(serIDs="character", con="TSodbcConnection"),
   definition= function(serIDs, con=getOption("TSconnection"),  
       vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
      TSdatesSQL(serIDs, con, vintage=vintage, panel=panel) )


setMethod("TSdescription",   signature(x="character", con="TSodbcConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        as.character(TSdescriptionSQL(x=x, con=con)) )

setMethod("TSdoc",   signature(x="character", con="TSodbcConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        as.character(TSdocSQL(x=x, con=con)) )

setMethod("TSlabel",   signature(x="character", con="TSodbcConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        as.character(TSlabelSQL(x=x, con=con)) )

setMethod("TSsource",   signature(x="character", con="TSodbcConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        TSsourceSQL(x=x, con=con))

setMethod("TSdelete", signature(serIDs="character", con="TSodbcConnection"),
     definition= function(serIDs, con=getOption("TSconnection"),  
     vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
   TSdeleteSQL(serIDs=serIDs, con=con, vintage=vintage, panel=panel) )

setMethod("TSvintages", 
   signature(con="TSodbcConnection"),
   definition=function(con) {
     if(!con@hasVintages) NULL else   
     sort(dbGetQuery(con,"SELECT  DISTINCT(vintage) FROM  vintages;" )$vintage)
     } )
