
TSdescriptionSQL <-  function(x=NULL, con=getOption("TSconnection"), 
       vintage=getOption("TSvintage"), panel=getOption("TSpanel"), 
       lang=getOption("TSlang")) {	    
            r <- dbGetQuery(con, paste("SELECT description", lang, 
	            "  FROM Meta ", setWhere(con, x, 
		        realVintage(con, vintage, x),
		        realPanel(con,panel)), ";", sep=""))[[1]]
	    # odbc converts NA to logical, but other db may return char
	    if(is.null(r) || is.na(r)|| ("NA" == r)) NA else r
	    }

TSdocSQL <-  function(x=NULL, con=getOption("TSconnection"), 
       vintage=getOption("TSvintage"), panel=getOption("TSpanel"), 
       lang=getOption("TSlang")) {
            if(1 < length(x)) stop("One series only for TSdoc")
            r <- dbGetQuery(con, paste("SELECT documentation", lang, 
	            "  FROM Meta ", setWhere(con, x, 
		        realVintage(con, vintage, x),
		        realPanel(con,panel)), ";", sep=""))[[1]]
	    # odbc converts NA to logical, but other db may return char
	    if(is.null(r) || is.na(r)|| ("NA" == r)) NA else r
	    }

TSlabelSQL <-  function(x=NULL, con=getOption("TSconnection"), 
       vintage=getOption("TSvintage"), panel=getOption("TSpanel"), 
       lang=getOption("TSlang")) {
            if(1 < length(x)) stop("One series only for TSlabel")
	    #  NOT YET
            #r <- dbGetQuery(con, paste("SELECT label", lang, 
	    #        "  FROM Meta ", setWhere(con, x, 
	    #            realVintage(con, vintage, x),
	    #	         realPanel(con,panel)), ";", sep=""))[[1]]
	    # odbc converts NA to logical, but other db may return char
	    #if(is.null(r) || is.na(r)|| ("NA" == r)) NA else r
	    NA
	    }

TSsourceSQL <-  function(x=NULL, con=getOption("TSconnection"), 
       vintage=getOption("TSvintage"), panel=getOption("TSpanel"), 
       lang=getOption("TSlang")) {
            if(1 < length(x)) stop("One series only for TSsource")
	    #  NOT YET
            #r <- dbGetQuery(con, paste("SELECT source", lang, 
	    #        "  FROM Meta ", setWhere(con, x, 
	    #            realVintage(con, vintage, x),
	    #	         realPanel(con,panel)), ";", sep=""))[[1]]
	    # odbc converts NA to logical, but other db may return char
	    #if(is.null(r) || is.na(r)|| ("NA" == r)) NA else r
	    NA
	    }

TSputSQL <- function(x, serIDs=seriesNames(x), con, Table=NULL,
       TSdescription.=TSdescription(x), TSdoc.=TSdoc(x),  TSlabel.=TSlabel(x),
       TSsource.=TSsource(x), 
       vintage=getOption("TSvintage"), panel=getOption("TSpanel")) {

  # so far I think this is generic to all SQL, but not extensively tested.
  # should rollback meta when data put fails
  #  (reversing order of M and P would mean data can change and 
  #    then meta fail, which seems worse.)
  panel <- realPanel(con,panel)
  # vintage <- realVintage(con,vintage, x) No. To put, vintage must already be real.
  # M does the write to Meta, P writes values to data table.

  # Column order could be specified after Meta, but this assumes  order
  #([vintage,] [ panel,] id, tbl, refperiod, description, documentation )
  M <- function(v, p, table, id, rP, desc, doc) {
    if(is.null(rP))   rP   <- "NA"
    if(is.null(desc)) desc <- "NA"
    if(is.null(doc))  doc  <- "NA"
    q <- "INSERT INTO Meta VALUES ('"
    if(!is.null(v))  q <- paste(q,v, "', '", sep="")
    if(!is.null(p))  q <- paste(q,p, "', '", sep="")
    q <- paste(q, id, "', '", table, "', '", 
                 rP, "', '", desc, "', '", doc, "') ;", sep="") 
    dbGetQuery(con, q) 
    # next if is RSQLite BUG work around
    if(is(con, "SQLiteConnection")) return(TRUE)
    0 == dbGetException(con)$errorNum
    }

  P <- function(p, tableV, columns, id, ...) {
    x <- list(...)
    vl <- id
    if(!is.null(p))  {
      vl <- paste(p,vl, sep=",")
      columns <- paste("panel,",columns, sep="")
      }
    for (i in 1:length(x)) vl <- paste(vl, x[[i]], sep="', '")
    # SQL wants NULL (not 'NULL') for NA
    vl <-   gsub("'NA'", "NULL", paste(vl,"'", sep=""))
    dbGetQuery(con, paste("DELETE FROM  ",tableV,
                     setWhere(con,id[1],NULL,panel), ";", sep="")) 
    for (i in seq(length(vl))){
      q <- paste("INSERT INTO ",tableV, " (", columns,
                     ") VALUES ('", vl[i], ") ", sep="")  
      dbGetQuery(con, q) 
      }
    # next if is RSQLite BUG work around
    if(is(con, "SQLiteConnection")) return(TRUE)
    0 == dbGetException(con)$errorNum
    }

  if (con@hasVintages) {
     hV <- TRUE
     vintage <- realVintage(con,vintage, serIDs)
     }
  else hV <- FALSE 

  # as.matrix(x) clobbers seriesNames(x)
  ids <-  serIDs 
  #ids <- gsub(" ", "", serIDs ) # remove spaces in id
  #if(! all( ids == serIDs)) warning("spaces removed from series names.")
  rP <- TSrefperiod(x)
  #N <- NROW(x) # could be periods() (or Tobs())
  # SHOULD PROBABLY HAVE METHODS FOR TS, ZOO, ITS, ... HERE
  if(is.ts(x)) {
    fr <-frequency(x) 
    y <- floor(time(x))
    p <- 1 + round(frequency(x) * (time(x) %% 1 ))
    x <- as.matrix(x)
    if(1==fr)   for (i in seq(nseries(x))) {  
	okM <-  M(vintage, panel, "A",  ids[i], rP[i], TSdescription.[i], TSdoc.[i])
        ok <- P(panel, tbNm(hV, "A", vintage), "id, year, v", ids[i], y, x[,i])
        }
    else if(2==fr)   for (i in seq(nseries(x))) {  
	okM <-  M(vintage, panel, "S",  ids[i], rP[i], TSdescription.[i], TSdoc.[i])
        ok <- P(panel, tbNm(hV, "S", vintage), "id, year, period, v", ids[i], y, p, x[,i])
        }
    else if(4==fr)   for (i in seq(nseries(x))) {  
	okM <-  M(vintage, panel, "Q",  ids[i], rP[i], TSdescription.[i], TSdoc.[i])
        ok <- P(panel, tbNm(hV, "Q", vintage), "id, year, period, v", ids[i], y, p, x[,i])
        }
    else if(12==fr)   for (i in seq(nseries(x))) {  
	okM <-  M(vintage, panel, "M",  ids[i], rP[i], TSdescription.[i], TSdoc.[i])
        ok <- P(panel, tbNm(hV, "M", vintage), "id, year, period, v", ids[i], y, p, x[,i])
        }
    else stop("ts frequency not supported.")  
    }
  else if (inherits(x, "zoo")) {
    ##require("zoo")
    #  might do better than this
    if (is.null(Table)) stop("Table must be specified for zoo series.")
    d <- time(x) 
    p <-  as.POSIXlt(d)
    x <- as.matrix(x)
    #  STILL NEED A S Q M  FOR ZOO
    if("W" == Table)   for (i in seq(nseries(x))) {  
	okM <-  M(vintage, panel, "W",  ids[i], rP[i], TSdescription.[i], TSdoc.[i])
        # period should be week of the year 1-52/3
	ok <- P(panel, tbNm(hV, "W", vintage), "id, date, period, v", ids[i], d, p$mday, x[,i])
        }
    else if("B" == Table)   for (i in seq(nseries(x))) {  
	okM <-  M(vintage, panel, "B",  ids[i], rP[i], TSdescription.[i], TSdoc.[i])
        # period should be business day of the year 1- ~260
	#  but need to map holidays
        ok <- P(panel, tbNm(hV, "B", vintage), "id, date, period, v", ids[i], d, p$yday, x[,i])
        }
    else if("D" == Table)   for (i in seq(nseries(x))) {  
	okM <-  M(vintage, panel, "D",  ids[i], rP[i], TSdescription.[i], TSdoc.[i])
        # period should be day of the year 1- 365/6
        ok <- P(panel, tbNm(hV, "D", vintage), "id, date, period, v", ids[i], d, p$yday, x[,i])
        }
    else if("U" == Table)   for (i in seq(nseries(x))) {  
	okM <-  M(vintage, panel, "U",  ids[i], rP[i], TSdescription.[i], TSdoc.[i])
        tz <- attr(p, "tzone")
	ok <- P(panel, tbNm(hV, "U", vintage), "id, date, period, v", ids[i], d, tz, p$mday,x[,i])
        }
    else if("T" == Table)   for (i in seq(nseries(x))) {  
	okM <-  M(vintage, panel, "T",  ids[i], rP[i], TSdescription.[i], TSdoc.[i])
        ok <- P(panel, tbNm(hV, "T", vintage), "id, date, v", ids[i], d, x[,i])
        }
    else if("I" == Table)   for (i in seq(nseries(x))) {  
	okM <-  M(vintage, panel, "I",  ids[i], rP[i], TSdescription.[i], TSdoc.[i])
        ok <- P(panel, tbNm(hV, "I", vintage), "id, date, v", ids[i], d, x[,i])
        }
    else stop("Table specification not supported.")  
     }
  else stop("Time series type not recognized.")  
  new("logicalId", ok & okM, 
        TSid=new("TSid", serIDs=serIDs, dbname=con@dbname, 
	         conType=class(con), hasVintages=con@hasVintages, hasPanels=con@hasPanels,
		 DateStamp=NA))
  }

TSdeleteSQL <- function(serIDs, con=getOption("TSconnection"),  
   vintage=getOption("TSvintage"), panel=getOption("TSpanel")) {
     for (i in seq(length(serIDs))) {
     	rv <- realVintage(con, vintage, i)
	rp <- realPanel(con,panel)
	where  <-  setWhere(con, serIDs[i], rv,   rp)
     	whereT <-  setWhere(con, serIDs[i], NULL, rp)
     	q <- dbGetQuery(con, paste("SELECT tbl  FROM Meta ",where, ";"))
	tbl <- tbNm(con@hasVintages, q$tbl, rv)
	if(0 != length(q)) {
     	   dbGetQuery(con, paste("DELETE FROM ", tbl, whereT, ";")) 
     	   dbGetQuery(con, paste("DELETE FROM Meta ",	where, ";")) 
     	   }
     	 }
     # could do better checking here
     new("logicalId", TRUE, 
          TSid=new("TSid", serIDs=serIDs, dbname=con@dbname, 
             conType=class(con), hasVintages=con@hasVintages, hasPanels=con@hasPanels, 
	     DateStamp=NA ))
     }

TSgetSQL <- function(serIDs, con, TSrepresentation=getOption("TSrepresentation"),
       tf=NULL, start=tfstart(tf), end=tfend(tf),
       names=NULL, 
       TSdescription=FALSE, TSdoc=FALSE, TSlabel=FALSE, TSsource=TRUE,
       vintage=getOption("TSvintage"), panel=getOption("TSpanel")) {
  # so far I think this is generic to all SQL.

  # default is ts if the table is in  c("A", "Q", "M","S") and zoo otherwise.
  # default is used for retrieval and any conversion is done after.

  if ( 1 < sum(c(length(serIDs), length(panel), length(vintage)) > 1))
   stop("Only one of serIDs, panel, or vintage can have length greater than 1.")

  if(is.null(names)) names <- 
         if ( length(panel)   > 1 )  panel   else
         if ( length(vintage) > 1 ) vintage  else  serIDs 

  panel <- realPanel(con,panel)
  # if vintage is a vector then serIDs needs to be expanded
  if ( 1 < length(vintage)) serIDs <- rep(serIDs, length(vintage))
  # next returns a vector of length equal serIDs
  if (con@hasVintages) {
     hV <- TRUE
     vintage <- realVintage(con,vintage, serIDs)
     }
  else hV <- FALSE 

  Q <- function(q) {# local function
      res <- dbGetQuery(con, q)
      if(any(dim(res) == 0)) stop("empty query result.")
      as.matrix(res)
      }

  mat <- desc <- doc <- label <- source <- rp <- NULL
  # if series are in "A", "Q", "M","S" use  ts otherwise zoo.
  for (i in seq(length(serIDs))) {
    where  <-  setWhere(con, serIDs[i], vintage[i], panel)
    whereT <-  setWhere(con, serIDs[i], NULL,       panel)
    for (j in seq(length(where))) {
    qq <- paste("SELECT tbl, refperiod  FROM Meta ",where[j], ";")
    q <- dbGetQuery(con, qq)
    if(0 == NROW(q$tbl))
       stop("Meta lookup for series ", serIDs[i], 
            " failed. (Result empty for query: ",
	     qq, ") Series does not exist on database.")

    if(1 < NROW(q$tbl)){
       warning("Meta lookup for series ", serIDs[i], 
            " returned multiple table entries ", q$tbl, ". query: ",
	     qq, ") Possible database corruption. Using the first table.")
	q <- q[1,]
	}

    if  (i == 1)  tbl <- q$tbl
    else if(q$tbl != tbl) 
       stop("Series must all have the same frequency or time representation.")

    rp <- c(rp, q$refperiod)

    if (tbl=="A") 
      {res <- Q(paste("SELECT year, v FROM ", tbNm(hV, "A", vintage[i]), 
                whereT[j], " order by year;"))
       r   <- ts(res[,2], start=c(res[1,1], 1), frequency=1) 
     }
    else if (tbl=="Q")  
      {res <- Q(paste("SELECT year, period, v FROM ", tbNm(hV, "Q", vintage[i]),
                whereT[j], " order by year, period;"))
       r   <- ts(res[,3], start=c(res[1,1:2]), frequency=4) 	 
      }
    else if (tbl=="M")
      {res <- Q(paste("SELECT year, period, v FROM ", tbNm(hV, "M", vintage[i]),
                whereT[j], " order by year, period;"))
       r   <- ts(res[,3], start=c(res[1,1:2]), frequency=12)	 
      }
    else if (tbl=="W") 
      {res <- Q(paste("SELECT date, period, v FROM ", tbNm(hV, "W", vintage[i]),
                whereT[j], " order by date;"))
       r   <- zoo::zoo(as.numeric(res[,3]), as.Date(res[,1]))
       # period is as.int(res[,2]) 	 
      }
    else if (tbl=="B") 
      {res <- Q(paste("SELECT date, period, v FROM ", tbNm(hV, "B", vintage[i]),
                whereT[j], " order by date;"))
       r   <- zoo::zoo(as.numeric(res[,3]), as.Date(res[,1]))
       # period is as.int(res[,2]) 	 
      }
    else if (tbl=="D")  
      {res <- Q(paste("SELECT date, period, v FROM ", tbNm(hV, "D", vintage[i]),
                whereT[j], " order by date;"))
       r   <- zoo::zoo(as.numeric(res[,3]), as.Date(res[,1]))
       # period is as.int(res[,2]) 	 
      }
    else if (tbl=="S")    
      {res <- Q(paste("SELECT year, period, v FROM ", tbNm(hV, "S", vintage[i]),
                whereT[j], " order by year, period;"))
       r   <- ts(res[,3], start=c(res[1,1:2]), frequency=2)	 
      }
    else if (tbl=="U")  
      {res <- Q(paste("SELECT date, tz, period, v FROM U ",tbNm(hV,"U",vintage[i]),
                whereT[j], " order by date;"))
       r   <- zoo::zoo(as.numeric(res[,4]), as.Date(res[,1]))
       # tz ? period is as.int(res[,3]) 	 
      }
    else if (tbl=="I")  
      {res <- Q(paste("SELECT date, v FROM ", tbNm(hV, "I", vintage[i]), 
                whereT[j], " order by date;"))
       r   <- zoo::zoo(as.numeric(res[,2]), as.Date(res[,1]))
      }
    else if (tbl=="T")  
      {res <- Q(paste("SELECT date, v FROM ", tbNm(hV, "T", vintage[i]), 
                whereT[j], " order by date;"))
       r   <- zoo::zoo(as.numeric(res[,2]), as.POSIXct(res[,1]))
      }
    else stop("Specified table not found.", 
              " (Internal TSdbi or database error likely.)",
	      " looking for series ", serIDs[i],
              " could not find tbl ", tbl)  
    
    if(TSdescription) desc <- c(desc,  TSdescription(serIDs[i],con) ) # where?
    if(TSdoc)         doc  <- c(doc,   TSdoc(serIDs[i],con) ) # where?
    if(TSlabel)       label<- c(label, TSlabel(serIDs[i],con) ) # where?
    if(TSsource)     source<- c(source,TSsource(serIDs[i],con) ) # where?

    mat <- tbind(mat, r)
    } # where[j]
    } # serID[i]

  mat <- tfwindow(mat, tf=tf, start=start, end=end)
  
  mat <- tframePlus::changeTSrepresentation(mat, TSrepresentation)

  if( (!all(is.na(rp))) && !all(rp == "	" ) ) TSrefperiod(mat) <- rp      

  seriesNames(mat) <- names

  TSmeta(mat) <- new("TSmeta", serIDs=serIDs, dbname=con@dbname, 
      conType=class(con), hasVintages=con@hasVintages, hasPanels=con@hasPanels, 
      DateStamp=NA, # bug in 2.7.0 =Sys.time(), 
      TSdescription=if(TSdescription) desc else  NA, 
      TSdoc        =if(TSdoc)          doc else  NA,
      TSlabel      =if(TSlabel)      label else  NA,
      TSsource     =if(TSsource)    source else  NA)
  mat
  }

TSdatesSQL <- function(serIDs, con,  
       vintage=getOption("TSvintage"), panel=getOption("TSpanel")) {
  # so far I think this is generic to all SQL, but untested.
  r  <- av <- tb <- rP <- NULL
  st <- en <- list()
  vintage <- realVintage(con, vintage, serIDs)
  panel   <- realPanel(con,panel)
  for (i in seq(length(serIDs))) {
    qq <- paste("SELECT id, tbl, refperiod  FROM Meta ", 
                    setWhere(con, serIDs[i], vintage[i], panel), ";", sep="")
    q <- dbGetQuery(con, qq)
    if(is.null(q) || NROW(q) == 0) {
        av <- c(av, FALSE)
	st <- append(st, list(NA))
	en <- append(en, list(NA))
	tb <- rbind(tb, NA)
	rP <- rbind(rP, NA)
	}
    else if(NROW(q) > 1){
      warning("More than one series with the same identifier. Possible database corruption.",
        " Meta lookup for series ", serIDs[i], " vintage ", vintage[i], 
            " returned multiple entries ", q, " for query: ",
	     qq, " Possible database corruption. Using the first table.")
	q <- q[1,]
	}
    else  {
      q2 <-  TSget(serIDs=serIDs[i], con, vintage=vintage, panel=panel)
      av <- c(av, TRUE)
      # paste(start(q2), collapse="-")
      st <- append(st, list(start(q2)))
      en <- append(en, list(end(q2)))
      tb <- rbind(tb, q$tb)
      rP <- rbind(rP, q$TSrefperiod)
     }
    }
  r <- serIDs
  attr(r, "TSdates") <- av
  attr(r, "start") <- st
  attr(r, "end")   <- en
  attr(r, "tbl")   <- tb
  attr(r, "TSrefperiod")   <- rP
  class(r) <- "TSdates"
  r
  }

#####################################
# these methods will generally not be needed by users, but is used in the test
# database setup. 

createTSdbTables <- function(con, index=FALSE){
 
 Texists <- function(a){
    if(dbExistsTable(con,a)) {
        warning("table ",a," exists. Not replacing it."); return(TRUE)}
    if(dbExistsTable(con,toupper(a))) {
        warning("table ",toupper(a)," exists. Not replacing it."); return(TRUE)}
    if(dbExistsTable(con,tolower(a))) {
        warning("table ",tolower(a)," exists. Not replacing it."); return(TRUE)}
    FALSE
    }
  
 if (!Texists("Meta")) {
  # Set up Metadata table  Meta.  
  dbGetQuery(con, "create table Meta (
     id 	 VARCHAR(40) NOT NULL,
     tbl	 CHAR(1), 
     refperiod   VARCHAR(10) default NULL,
     description   TEXT,
     documentation     TEXT
     );")
  dbGetQuery(con, "CREATE UNIQUE INDEX Metaindex_id ON Meta (id);")
  dbGetQuery(con, "CREATE INDEX Metaindex_tbl ON Meta (tbl);")
  }
  
 if (!Texists("A")) {
  # Set up annual table  A.
  dbGetQuery(con, "create table A (
     id 	VARCHAR(40),
     year	INT,
     v    double precision DEFAULT  NULL
     );")
  if (index) dbGetQuery(con, "CREATE INDEX Aindex_id     ON A (id);")
  if (index) dbGetQuery(con, "CREATE INDEX Aindex_year   ON A (year);")
  }
  
 if (!Texists("B")) {
  # Set up business table  B . 
  dbGetQuery(con, "create table B (
     id 	VARCHAR(40),
     date	DATE,
     period	INT,
     v    double precision DEFAULT  NULL
     );")
  if (index) dbGetQuery(con, "CREATE INDEX Bindex_id     ON B (id);")
  if (index) dbGetQuery(con, "CREATE INDEX Bindex_date   ON B (date);")
  if (index) dbGetQuery(con, "CREATE INDEX Bindex_period ON B (period);")
  }
  
 if (!Texists("D")) {
  # Set up daily table  D .  
  dbGetQuery(con, "create table D (
     id 	VARCHAR(40),
     date	DATE,
     period	INT,
     v    double precision DEFAULT  NULL
     );")
  if (index) dbGetQuery(con, "CREATE INDEX Dindex_id     ON D (id);")
  if (index) dbGetQuery(con, "CREATE INDEX Dindex_date   ON D (date);")
  if (index) dbGetQuery(con, "CREATE INDEX Dindex_period ON D (period);")
  }
  
 if (!Texists("M")) {
  # Set up monthly table  M.
  dbGetQuery(con, "create table M (
     id 	VARCHAR(40),
     year	INT,
     period	INT,  
     v    double precision DEFAULT  NULL
     );")
  if (index) dbGetQuery(con, "CREATE INDEX Mindex_id     ON M (id);")
  if (index) dbGetQuery(con, "CREATE INDEX Mindex_year   ON M (year);")
  if (index) dbGetQuery(con, "CREATE INDEX Mindex_period ON M (period);")
  }
  
 if (!Texists("U")) {
  # Set up minutely table  U .
  # tz not tested yet. Not sure about period.  
  dbGetQuery(con, "create table U (
     id 	VARCHAR(40),
     date	TIMESTAMP,
     tz 	VARCHAR(4),    
     period	INT,	     
     v    double precision DEFAULT  NULL
     );")
  if (index) dbGetQuery(con, "CREATE INDEX Uindex_id     ON U (id);")
  if (index) dbGetQuery(con, "CREATE INDEX Uindex_date   ON U (date);")
  if (index) dbGetQuery(con, "CREATE INDEX Uindex_period ON U (period);")
  }
  
 if (!Texists("Q")) {
  # Set up quarterly table  Q.
  dbGetQuery(con, "create table Q (
     id 	VARCHAR(40),
     year	INT,
     period	INT,  
     v    double precision DEFAULT  NULL
     );")
  if (index) dbGetQuery(con, "CREATE INDEX Qindex_id     ON Q (id);")
  if (index) dbGetQuery(con, "CREATE INDEX Qindex_year   ON Q (year);")
  if (index) dbGetQuery(con, "CREATE INDEX Qindex_period ON Q (period);")
  }
  
 if (!Texists("S")) {
  # Set up semiannual table  S .
  dbGetQuery(con, "create table S (
     id 	VARCHAR(40),
     year	INT,
     period	INT,  
     v    double precision DEFAULT  NULL
     );")
  if (index) dbGetQuery(con, "CREATE INDEX Sindex_id     ON S (id);")
  if (index) dbGetQuery(con, "CREATE INDEX Sindex_year   ON S (year);")
  if (index) dbGetQuery(con, "CREATE INDEX Sindex_period ON S (period);")
  }
  
 if (!Texists("W")) {
  # Set up weekly table  W .
  dbGetQuery(con, "create table W (
     id 	VARCHAR(40),
     date	DATE,
     period	INT,  
     v    double precision DEFAULT  NULL
     );")
  if (index) dbGetQuery(con, "CREATE INDEX Windex_id     ON W (id);")
  if (index) dbGetQuery(con, "CREATE INDEX Windex_date   ON W (date);")
  if (index) dbGetQuery(con, "CREATE INDEX Windex_period ON W (period);")
  }
  
 if (!Texists("I")) {
  # Set up irregular date table  I .
  dbGetQuery(con, "create table I (
     id 	VARCHAR(40),
     date	DATE,
   v    double precision DEFAULT  NULL
   );")
  if (index) dbGetQuery(con, "CREATE INDEX Iindex_id     ON I (id);")
  if (index) dbGetQuery(con, "CREATE INDEX Iindex_date   ON I (date);")
  }
  
 if (!Texists("T")) {
  # Set up irregular date-time table  T .
  # IT WOULD BE NICE TO HAVE TZ HERE TOO
  dbGetQuery(con, "create table T (
     id    VARCHAR(40),
     date  TIMESTAMP,
     v     double precision DEFAULT  NULL
     );")
  if (index) dbGetQuery(con, "CREATE INDEX Tindex_id   ON T (id);")
  if (index) dbGetQuery(con, "CREATE INDEX Tindex_date ON T (date);")
  }
    
  #  This is generic sql way to get table info. (eg. table A )
  #   but it requires read privileges on INFORMATION_SCHEMA.Columns
  #   which the user may not have.
  #
  #if( "SQLiteConnection" != class(con)) { #SQLite does not seem to support this
  #  z <- try(
  #  dbGetQuery(con, "SELECT COLUMN_NAME, COLUMN_DEFAULT, COLLATION_NAME, DATA_TYPE,
  #	   CHARACTER_SET_NAME, CHARACTER_MAXIMUM_LENGTH, NUMERIC_PRECISION
  #	      FROM INFORMATION_SCHEMA.Columns WHERE table_name='A' ;"), silent=TRUE )
   # if(inherits(z, "try-error")) cat(
  #    "INFORMATION_SCHEMA query problem. (User may not have permission.)\n")
  #  else print(z)
  # }

  cat("   tables:\n")
  dbListTables(con)
  }
   
removeTSdbTables <- function(con, yesIknowWhatIamDoing=FALSE, ToLower=FALSE){

  if(!yesIknowWhatIamDoing) 
       stop("You need to know what you are doing before using this function!!")
   
  # "dropTStable" works around the
  # problem that different db engines treat capitalized table names differently.
  # e.g. MySQL uses table name Meta while Posgresql converts to meta.

  dropTStable <- function(Table){
    if (ToLower) Table <- tolower(Table)
    if(dbExistsTable(con, Table)) dbRemoveTable(con, Table)
    }

  dropTStable("Meta")   
  dropTStable("A")
  dropTStable("B")
  dropTStable("D")
  dropTStable("M")
  dropTStable("U")
  dropTStable("Q")
  dropTStable("S")
  dropTStable("W")
  dropTStable("I")
  dropTStable("T")
  invisible(TRUE)
  }
  
   
########## internal utilities to construct WHERE and table name#########

tbNm <- function(hasVintages, tbl, rV) {
	if(hasVintages) tbl <- paste(tbl, rV, sep="_")
	# map - to _ in table names
	gsub("-", "_",tbl)
	}

realVintage <- function(con, vintage, serIDs) {
   # replace alias with canonical name if necessary
   # assuming if vintage is a vector then serIDs has already been expanded to
   # the same length. Otherwise, vintage should be null or a scalar which is
   # expanded to the length of serIDs. 
   if(!con@hasVintages) return(vintage) #usually NULL in this case 
   
   if(is.null(vintage)) vintage <- "current"
   else if( 1 == length(vintage)) vintage <- rep(vintage, length(serIDs))
   
   rV <- NULL
   for (i in seq(length(serIDs))){
     q <- paste("SELECT vintage  FROM vintageAlias WHERE alias ='",
              vintage[i],"' AND id = '",serIDs[i],"';", sep="") 
     r <- dbGetQuery(con,q )$vintage
     # if alias result is empty check null id which applies to all series.
     if (0== NROW(r)) {
        q <- paste( "SELECT vintage  FROM vintageAlias WHERE alias='", 
	            vintage[i],"' ;", sep="")
        r <- dbGetQuery(con,q )$vintage
        }
     # if alias result is still empty assume vintage is already the real one.
     if (0== NROW(r)) r <- vintage[i]
     rV <- c(rV, r)
     }
   rV
   }

realPanel <- function(con, panel) {
   # replace alias with canonical name if necessary
   if(!con@hasPanels) return(panel) #usually NULL in this case
   if(is.null(panel)) stop("panel must be specified")
   for (i in seq(length(panel))){
     q <- paste("SELECT panel  FROM panelAlias WHERE alias='",
                panel[i],"';", sep="") 

     r <- dbGetQuery(con,q )$panel
     # if alias result is empty assume panel is already the real one.
     if (0== NROW(r)) r <- panel[i]
     rP <- c(rP, r)
     }
   rP
   }

setWhere <- function(con, serIDs, realVintage, realPanel) {
   # serIDs must be a scale for this function
   # Calls for Meta will pass realVintage, but data tables will set 
   #  realVintage=NULL and append vintage to table name.
   where <-  paste(" WHERE id = '", serIDs, "'", sep="")
   if(!is.null(realVintage))
      where <- paste(where, " AND vintage='", realVintage, "'", sep="")
   if(!is.null(realPanel)) 
      where <- paste(where, " AND panel='",   realPanel,   "'", sep="")
   where
   }

