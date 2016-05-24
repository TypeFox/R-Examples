SQLite <- function(...){
  args <- pmatch(names(list(...)), names(formals(RSQLite::SQLite)) )
  do.call(RSQLite::SQLite, list(...)[!is.na(args)])
  }

setClass("TSSQLiteConnection", contains=c("SQLiteConnection","conType", "TSdb")) 

#setAs("TSSQLiteConnection", "integer", 
#  def=getMethod("coerce", c("dbObjectId","integer"))) 

setGeneric("print")

setMethod("print", "TSSQLiteConnection", function(x, ...) {
    print(x@TSdb)
    })

#class(getExportedValue("RSQLite", "SQLite")()) is "SQLiteDriver"

setMethod("TSconnect",   signature(q="SQLiteConnection", dbname="missing"),
   definition=function(q, dbname, ...) {
        con <- q
	if(!dbExistsTable(conn=con, "Meta"))
	  stop("The database does not appear to be a TS database,")
	nm <- dbGetQuery(con, "PRAGMA database_list;")
	nm <- nm$file["main" == nm$name]
	new("TSSQLiteConnection" , con, dbname=nm, 
	       hasVintages=dbExistsTable(conn=con, "vintageAlias"), 
	       hasPanels  =dbExistsTable(conn=con, "panels"))
	})

setMethod("TSput", signature(x="ANY", serIDs="character", con="TSSQLiteConnection"),
   definition= function(x, serIDs=seriesNames(x), con=getOption("TSconnection"), Table=NULL, 
       TSdescription.=TSdescription(x), TSdoc.=TSdoc(x),   TSlabel.=TSlabel(x),
         TSsource.=TSsource(x),
       vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...) 
  TSputSQL(x, serIDs, con, Table=Table, 
     TSdescription.=TSdescription., TSdoc.=TSdoc., TSlabel.=TSlabel.,
     TSsource.=TSsource.,
  vintage=vintage, panel=panel) )

setMethod("TSget", signature(serIDs="character", con="TSSQLiteConnection"),
   definition= function(serIDs, con=getOption("TSconnection"), 
       TSrepresentation=options()$TSrepresentation,
       tf=NULL, start=tfstart(tf), end=tfend(tf), names=NULL, 
       TSdescription=FALSE, TSdoc=FALSE, TSlabel=FALSE, TSsource=TRUE,
       vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
    TSgetSQL(serIDs, con, TSrepresentation=TSrepresentation,
       tf=tf, start=start, end=end, names=names, 
       TSdescription=TSdescription, TSdoc=TSdoc, TSlabel=TSlabel,
         TSsource=TSsource,
       vintage=vintage, panel=panel) )

setMethod("TSdates", signature(serIDs="character", con="TSSQLiteConnection"),
   definition= function(serIDs, con=getOption("TSconnection"),  
       vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
     TSdatesSQL(serIDs, con, vintage=vintage, panel=panel) )


setMethod("TSdescription",   signature(x="character", con="TSSQLiteConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        TSdescriptionSQL(x=x, con=con) )

setMethod("TSdoc",   signature(x="character", con="TSSQLiteConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        TSdocSQL(x=x, con=con) )

setMethod("TSlabel",   signature(x="character", con="TSSQLiteConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        TSlabelSQL(x=x, con=con) )

setMethod("TSsource",   signature(x="character", con="TSSQLiteConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        TSsourceSQL(x=x, con=con) )

setMethod("TSdelete", signature(serIDs="character", con="TSSQLiteConnection"),
     definition= function(serIDs, con=getOption("TSconnection"),  
     vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
  TSdeleteSQL(serIDs=serIDs, con=con, vintage=vintage, panel=panel) )


setMethod("TSvintages", 
   signature(con="TSSQLiteConnection"),
   definition=function(con) {
     if(!con@hasVintages) NULL else   
     sort(dbGetQuery(con,"SELECT  DISTINCT(vintage) FROM  vintages;" )$vintage)
     } )
