PostgreSQL <- function(...){
  args <- pmatch(names(list(...)), names(formals(RPostgreSQL::PostgreSQL)) )
  do.call(RPostgreSQL::PostgreSQL, list(...)[!is.na(args)])
  }

setClass("TSPostgreSQLConnection", 
   contains=c("PostgreSQLConnection", "conType", "TSdb")) 



#setMethod("print", "TSPostgreSQLConnection", function(x, ...) {
#   print(x@TSdb)
#   })

#class(getExportedValue("RPostgreSQL", "PostgreSQL")()) is "RPostgreSQL"

setMethod("TSconnect",   signature(q="PostgreSQLConnection", dbname="missing"),
   definition=function(q, dbname, host=
    if(!is.null(Sys.getenv("PGHOST"))) Sys.getenv("PGHOST") else "localhost", ...) {
        con <- q
	nm <- as.character(dbGetQuery(conn=con, "SELECT  current_database();"))
	if(0 == length(dbListTables(con))){
	  dbDisconnect(con)
          stop("Database ",dbname," has no tables.")
	  }
	if(!dbExistsTable(con, "meta")){ # meta not Meta for PostgreSQL
	  dbDisconnect(con)
          stop("Database ",dbname," does not appear to be a TS database.")
	  }
  	new("TSPostgreSQLConnection" , con, dbname=nm, 
 	       hasVintages=dbExistsTable(con, "vintagealias"),  # vintagealias not vintageAlias for PostgreSQL
 	       hasPanels  =dbExistsTable(con, "panels")) 
	})

setMethod("TSput",   signature(x="ANY", serIDs="character", con="TSPostgreSQLConnection"),
   definition= function(x, serIDs, con=getOption("TSconnection"), Table=NULL, 
       TSdescription.=TSdescription(x), TSdoc.=TSdoc(x), TSlabel.=TSlabel(x),
       TSsource.=TSsource(x),  
       vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
   TSputSQL(x, serIDs, con, Table=Table, 
   TSdescription.=TSdescription., TSdoc.=TSdoc., TSlabel.=TSlabel.,
    TSsource.=TSsource.,
   vintage=vintage, panel=panel) )

setMethod("TSget",   signature(serIDs="character", con="TSPostgreSQLConnection"),
   definition= function(serIDs, con=getOption("TSconnection"), 
       TSrepresentation=getOption("TSrepresentation"),
       tf=NULL, start=tfstart(tf), end=tfend(tf), names=NULL, 
       TSdescription=FALSE, TSdoc=FALSE, TSlabel=FALSE, TSsource=TRUE,
       vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
     TSgetSQL(serIDs, con, TSrepresentation=TSrepresentation,
       tf=tf, start=start, end=end, names=names, 
       TSdescription=TSdescription, TSdoc=TSdoc, TSlabel=TSlabel,
         TSsource=TSsource,
       vintage=vintage, panel=panel) )

setMethod("TSdates",    signature(serIDs="character", con="TSPostgreSQLConnection"),
   definition= function(serIDs, con=getOption("TSconnection"),  
       vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
      TSdatesSQL(serIDs, con, vintage=vintage, panel=panel) )


setMethod("TSdescription",   signature(x="character", con="TSPostgreSQLConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        TSdescriptionSQL(x=x, con=con) )

setMethod("TSdoc",   signature(x="character", con="TSPostgreSQLConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        TSdocSQL(x=x, con=con) )

setMethod("TSlabel",   signature(x="character", con="TSPostgreSQLConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        TSlabelSQL(x=x, con=con) )

setMethod("TSsource",   signature(x="character", con="TSPostgreSQLConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        TSsourceSQL(x=x, con=con) )

setMethod("TSdelete", signature(serIDs="character", con="TSPostgreSQLConnection"),
   definition= function(serIDs, con=getOption("TSconnection"),  
   vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
       TSdeleteSQL(serIDs=serIDs, con=con, vintage=vintage, panel=panel) )

setMethod("TSvintages", 
   signature(con="TSPostgreSQLConnection"),
   definition=function(con) {
     if(!con@hasVintages) NULL else   
     sort(dbGetQuery(con,"SELECT  DISTINCT(vintage) FROM  vintages;" )$vintage)
     } )
