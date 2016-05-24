MySQL <- function(...){
  args <- pmatch(names(list(...)), names(formals(RMySQL::MySQL)) )
  do.call(RMySQL::MySQL, list(...)[!is.na(args)])
  }

setClass("TSMySQLConnection", contains=c("MySQLConnection", "conType", "TSdb"))

# works  with both pre-release 2015.1 and older versions of RMySQL
# but setAs really should not be needed
#setAs("TSMySQLConnection", "integer",
# def = function(from) as(slot(from,"Id"), "integer")
# )

# in which case we need 
#new("TSMySQLConnection" , con, drv="MySQL", dbname=dbname, 
#  	       hasVintages=dbExistsTable(con, "vintages"), 
#  	       hasPanels  =dbExistsTable(con, "panels")) 

#setMethod("print", "TSMySQLConnection", function(x, ...) {
#    print(x@TSdb)
#    })

#class(getExportedValue("RMySQL", "MySQL")()) is "MySQLDriver"

setMethod("TSconnect",   signature(q="MySQLConnection", dbname="missing"),
   definition=function(q, dbname, ...) {
        con <- q
	nm <- as.character(dbGetQuery(conn=con, "SELECT DATABASE();"))
	if(0 == length(dbListTables(con))){
	  dbDisconnect(con)
          stop("Database ",nm," has no tables.")
	  }
	if(!dbExistsTable(con, "Meta")){
	  dbDisconnect(con)
          stop("Database ",nm," does not appear to be a TS database.")
	  }
	new("TSMySQLConnection" , con, dbname=nm, 
  	       hasVintages=dbExistsTable(con, "vintageAlias"), 
  	       hasPanels  =dbExistsTable(con, "panels")) 
	})

setMethod("TSput",   signature(x="ANY", serIDs="character", con="TSMySQLConnection"),
   definition= function(x, serIDs, con=getOption("TSconnection"), Table=NULL, 
       TSdescription.=TSdescription(x), TSdoc.=TSdoc(x), TSlabel.=TSlabel(x),
       TSsource.=TSsource(x),  
       vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
       TSputSQL(x, serIDs, con, Table=Table, 
       TSdescription.=TSdescription., TSdoc.=TSdoc., TSlabel.=TSlabel.,
       TSsource.=TSsource.,
       vintage=vintage, panel=panel) )

setMethod("TSget",   signature(serIDs="character", con="TSMySQLConnection"),
   definition= function(serIDs, con=getOption("TSconnection"), 
       TSrepresentation=options()$TSrepresentation,
       tf=NULL, start=tfstart(tf), end=tfend(tf),
       names=NULL, TSdescription=FALSE, TSdoc=FALSE,TSlabel=FALSE,TSsource=TRUE,
       vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
       TSgetSQL(serIDs, con, TSrepresentation=TSrepresentation,
       tf=tf, start=start, end=end,
       names=names, TSdescription=TSdescription, TSdoc=TSdoc, TSlabel=TSlabel,
       TSsource=TSsource,
       vintage=vintage, panel=panel) )

setMethod("TSdates",    signature(serIDs="character", con="TSMySQLConnection"),
   definition= function(serIDs, con=getOption("TSconnection"),  
       vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
       TSdatesSQL(serIDs, con, vintage=vintage, panel=panel) )


setMethod("TSdescription",   signature(x="character", con="TSMySQLConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
         TSdescriptionSQL(x=x, con=con) )

setMethod("TSdoc",   signature(x="character", con="TSMySQLConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        TSdocSQL(x=x, con=con) )

setMethod("TSlabel",   signature(x="character", con="TSMySQLConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        TSlabelSQL(x=x, con=con) )

setMethod("TSsource",   signature(x="character", con="TSMySQLConnection"),
   definition= function(x, con=getOption("TSconnection"), ...)
        TSsourceSQL(x=x, con=con) )

setMethod("TSdelete", 
   signature(serIDs="character", con="TSMySQLConnection", vintage="ANY", panel="ANY"),
   definition= function(serIDs, con=getOption("TSconnection"),  
   vintage=getOption("TSvintage"), panel=getOption("TSpanel"), ...)
   TSdeleteSQL(serIDs=serIDs, con=con, vintage=vintage, panel=panel) )

setMethod("TSvintages", 
   signature(con="TSMySQLConnection"),
   definition=function(con) {
     if(!con@hasVintages) NULL else   
     sort(dbGetQuery(con,"SELECT  DISTINCT(vintage) FROM  vintages;" )$vintage)
     } )
