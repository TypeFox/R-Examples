###############################################################################
## Project: ETLUtils
## Content: Utilities that work for ff objects
##
## Author: jan
## Creation date: Mar 25, 2012, 10:06:53 PM
## File: ffsql.R
###############################################################################
	


#' Read data from a DBI connection into an ffdf.
#' 
#' Read data from a DBI connection into an \code{\link[ff]{ffdf}}. This can for example be used to import
#' large datasets from Oracle, SQLite, MySQL, PostgreSQL, Hive or other SQL databases into R. \cr
#' 
#' Opens up the DBI connection using \code{DBI::dbConnect}, sends the query using \code{DBI::dbSendQuery} and \code{DBI::fetch}-es 
#' the results in batches of next.rows rows. Heavily borrowed from \code{\link[ff]{read.table.ffdf}} 
#'
#' @param query the SQL query to execute on the DBI connection
#' @param dbConnect.args a list of arguments to pass to DBI's \code{\link[DBI]{dbConnect}} (like drv, dbname, username, password). See the examples.
#' @param dbSendQuery.args a list containing database-specific parameters which will be passed to to pass to \code{\link[DBI]{dbSendQuery}}.
#' Defaults to an empty list.  
#' @param fetch.args a list containing optional database-specific parameters which will be passed to to pass to \code{\link[DBI]{fetch}}. 
#' Defaults to an empty list.          
#' @param x NULL or an optional ffdf object to which the read records are appended. 
#' See documentation in read.table.ffdf for more details and the example below.
#' @param nrows Number of rows to read from the query resultset. Default value of -1 reads in all rows.
#' @param first.rows chunk size (rows) to read for first chunk from the query resultset
#' @param next.rows chunk size (rows) to read sequentially for subsequent chunks from the query resultset. Currently, this must be specified.
#' @param levels optional specification of factor levels. A list with as names the names the columns of the data.frame 
#' fetched in the first.rows, containing levels of the factors. 
#' @param appendLevels logical. A vector of permissions to expand levels for factor columns. See documentation in \code{\link[ff]{read.table.ffdf}} for more details.
#' @param asffdf_args further arguments passed to \code{\link[ff]{as.ffdf}} (ignored if 'x' gives an ffdf object )
#' @param BATCHBYTES integer: bytes allowed for the size of the data.frame storing the result of reading one chunk. 
#' See documentation in \code{\link[ff]{read.table.ffdf}} for more details.
#' @param VERBOSE logical: TRUE to verbose timings for each processed chunk (default FALSE).
#' @param colClasses See documentation in \code{\link[ff]{read.table.ffdf}}
#' @param transFUN function applied to the data frame after each chunk is retreived by \code{\link[DBI]{fetch}}
#' @param ... optional parameters passed on to transFUN
#' @return 
#' An ffdf object unless the query returns zero records in which case the function will return the data.frame
#' returned by \code{\link[DBI]{fetch}} and possibly transFUN. 
#' @export
#' @seealso \code{\link[ff]{read.table.ffdf}, \link{read.odbc.ffdf}}
#' @examples
#' require(ff)
#' 
#' ##
#' ## Example query using data in sqlite
#' ##
#' require(RSQLite)
#' dbfile <- system.file("smalldb.sqlite3", package="ETLUtils")
#' drv <- dbDriver("SQLite")
#' query <- "select * from testdata limit 10000"
#' x <- read.dbi.ffdf(query = query, dbConnect.args = list(drv = drv, dbname = dbfile), 
#' first.rows = 100, next.rows = 1000, VERBOSE=TRUE)
#' class(x)
#' x[1:10, ]
#' 
#' ## show it is the same as getting the data directly using RSQLite 
#' ## apart from characters which are factors in ffdf objects
#' directly <- dbGetQuery(dbConnect(drv = drv, dbname = dbfile), query)
#' directly <- as.data.frame(as.list(directly), stringsAsFactors=TRUE)
#' all.equal(x[,], directly)
#'
#' ## show how to use the transFUN argument to transform the data before saving into the ffdf
#' ## and shows the use of the levels argument
#' query <- "select * from testdata limit 10"
#' x <- read.dbi.ffdf(query = query, dbConnect.args = list(drv = drv, dbname = dbfile), 
#' first.rows = 100, next.rows = 1000, VERBOSE=TRUE, levels = list(a = rev(LETTERS)),
#' transFUN = function(x, subtractdays){
#' 	x$b <- as.Date(x$b)
#' 	x$b.subtractdaysago <- x$b - subtractdays
#' 	x
#' }, subtractdays=7)
#' class(x)
#' x[1:10, ]
#' ## remark that the levels of column a are reversed due to specifying the levels argument correctly
#' levels(x$a)
#' 
#' ## show how to append data to an existing ffdf object 
#' transformexample <- function(x, subtractdays){
#' 	x$b <- as.Date(x$b)
#' 	x$b.subtractdaysago <- x$b - subtractdays
#' 	x
#' }
#' dim(x)
#' x[,]
#' combined <- read.dbi.ffdf(query = query, 
#'  dbConnect.args = list(drv = drv, dbname = dbfile), 
#'  first.rows = 100, next.rows = 1000, x = x, VERBOSE=TRUE, 
#'  transFUN = transformexample, subtractdays=1000)
#' dim(combined)
#' combined[,]
#' 
#' ##
#' ## Example query using ROracle. Do try this at home with some larger data :)
#' ##
#' \dontrun{
#' require(ROracle)
#' query <- "select OWNER, TABLE_NAME, TABLESPACE_NAME, NUM_ROWS, LAST_ANALYZED from all_all_tables" 
#' x <- read.dbi.ffdf(query=query,
#' dbConnect.args = list(drv = dbDriver("Oracle"), 
#' user = "YourUser", password = "YourPassword", dbname = "Mydatabase"),
#' first.rows = 100, next.rows = 50000, nrows = -1, VERBOSE=TRUE)
#' }
read.dbi.ffdf <- function(
		query = NULL,
		dbConnect.args = list(drv=NULL, dbname = NULL, username = "", password = ""), 
		dbSendQuery.args = list(), 
		fetch.args = list(), 		
		x = NULL, nrows = -1, 
		first.rows = NULL, next.rows = NULL, levels = NULL, appendLevels = TRUE, 
		asffdf_args = list(), BATCHBYTES = getOption("ffbatchbytes"), VERBOSE = FALSE, colClasses = NULL, 
		transFUN = NULL, ...){
	cleanupConnection <- function(x){
		if("resultset" %in% names(x)){
			DBI::dbClearResult(res=x$resultset)
		}
		if("channel" %in% names(x)){
		  DBI::dbDisconnect(x$channel)
		}
	}
	dbiinfo <- list()
	
	##
	## Connect to database
	##
	dbiinfo$channel <- do.call('dbConnect', dbConnect.args)
	on.exit(cleanupConnection(dbiinfo))
	
	append <- !is.null(x)
	if(append && !inherits(x, "ffdf")){
		stop("only ffdf objects can be used for appending (and skipping the first.row chunk)")
	}
	if(VERBOSE){
		read.time <- 0
		write.time <- 0
	}
	
	nrows <- as.integer(nrows)
	N <- 0L	
	if(!append){		
		if(VERBOSE){
			cat("read.dbi.ffdf ", N+1L, "..", sep="")
			read.start <- proc.time()[3]
		}		
		if(is.null(colClasses)){
			colClasses <- NA
		}		
		## Specify the size of the first number of rows to import
		if (is.null(first.rows)){
			if (is.null(next.rows)){
				first.rows <- 1e3L
			}else{
				first.rows <- as.integer(next.rows)
			}				
		}else{
			first.rows <- as.integer(first.rows)
		}		
		if (nrows>=0L && nrows<first.rows){
			first.rows <- nrows
		}		
		if (first.rows==1){
			stop("first.rows must not be 1")
		}
		##
		## Set up arguments for odbcQuery and sqlGetResults functions
		##
		dbSendQuery.args$conn <- dbiinfo$channel
		dbSendQuery.args$statement <- query
		##
		## Launch the query
		##
		dbiinfo$resultset <- do.call('dbSendQuery', dbSendQuery.args)
				
		##
		## Read in the first chunk
		##
		fetch.args$res <- dbiinfo$resultset
		fetch.args$n <- first.rows			
		if(is.null(transFUN)){
			dat <- convertCharacterToFactor(do.call("fetch", fetch.args))
		}else{
			dat <- convertCharacterToFactor(transFUN(do.call("fetch", fetch.args), ...))
		}				
		n <- nrow(dat)
		N <- n
		## Recode levels if levels are specified
		if(!is.null(levels)){
			cnam <- colnames(dat)
			lnam <- names(levels)
			if(!is.list(levels) || is.null(lnam) || any(is.na(match(lnam,cnam)))){
				stop("levels must be a list with names matching column names of the first data.frame read")
			}				
			for (i in lnam){
				dat[[i]] <- recodeLevels(dat[[i]], levels[[i]])
			}
		}		
		if (VERBOSE){
			write.start <- proc.time()[3]
			read.time <- read.time + (write.start - read.start)
			cat(N, " (", n, ")  dbi-read=", round(write.start-read.start, 3), "sec", sep="")
		}
		##
		## If there are zero rows in the query, return a data.frame instead of an ffdf object
		##
		if(nrow(dat) == 0){
			if (VERBOSE){
				cat(" query returned 0 records\n", sep="")
				cat(" dbi-read=", round(read.time, 3), "sec  TOTAL=", round(read.time, 3), "sec\n", sep="")
			}				
			return(dat)
		}
		
		##
		## If there are non-zero rows in the query, convert to a filebased ffdf object
		##
		x <- do.call("as.ffdf", c(list(dat), asffdf_args))		
		##
		## Fix ordered factors
		##
		colClasses <- repnam(colClasses, colnames(x), default=NA)
		i.fix <- seq.int(length.out=ncol(dat))[!is.na(match(colClasses, "ordered"))]
		for (i in i.fix){
			virtual(x[[i]])$ramclass <- c("ordered","factor")
		}		
		if (VERBOSE){
			write.stop <- proc.time()[3]
			write.time <- write.time + (write.stop - write.start)
			cat(" ffdf-write=", round(write.stop-write.start, 3), "sec\n", sep="")
		}
	}
	##
	## Fetch the other rows
	##
	if (append || N==first.rows){
		##
		## Define the number of rows to fetch in each run depending on the BATCHBYTES maximum
		##
		k <- ncol(x)
		col.names <- colnames(x)
		colClasses <- sapply(seq.int(length.out=ncol(x)), function(i)colClass(x[[i]]))		
		if (is.null(next.rows)){
			recordsize <- sum(.rambytes[vmode(x)])
			next.rows <- BATCHBYTES %/% recordsize
			if (next.rows<1L){
				next.rows <- 1L
				warning("single record does not fit into BATCHBYTES")
			}
		}else{
			next.rows <- as.integer(next.rows)
		}		
		next.nrows <- next.rows
		
		appendLevels <- repnam(appendLevels, col.names, default=TRUE)
		if(any(appendLevels)){
			i.fac <- seq.int(length.out=k)
			if(append == TRUE){
				i.fac <- i.fac[appendLevels & sapply(i.fac, function(i) is.factor(x[[i]]))]
			}else{
				i.fac <- i.fac[appendLevels & sapply(i.fac, function(i) is.factor(dat[[i]]))]	
			}
			
		}	
		if(append == TRUE){
			## Set up arguments for odbcQuery and sqlGetResults functions & Launch the query
			dbSendQuery.args$conn <- dbiinfo$channel
			dbSendQuery.args$statement <- query			
			dbiinfo$resultset <- do.call('dbSendQuery', dbSendQuery.args)
		}
		##
		## Loop - fetching the subsequent data in batches
		##
		while(TRUE){
			if (nrows>=0L && N+next.rows > nrows){
				next.nrows <- nrows - N
			}				
			if (next.nrows<1L){
				break
			}			
			if (VERBOSE){
				cat("read.dbi.ffdf ", N+1L,"..", sep="")
				read.start <- proc.time()[3]
			}
			##
			## read in subsequent chunk
			##
			fetch.args$res <- dbiinfo$resultset
			fetch.args$n <- next.nrows
			if (is.null(transFUN)){
				dat <- convertCharacterToFactor(do.call("fetch", fetch.args))
			}else{
				dat <- convertCharacterToFactor(transFUN(do.call("fetch", fetch.args), ...))
			}			
			n <- nrow(dat)
			N <- N + n
			if (VERBOSE){
				write.start <- proc.time()[3]
				read.time <- read.time + (write.start - read.start)
				cat(N, " (", n, ")  dbi-read=", round(write.start-read.start, 3), "sec", sep="")
			}
			## Nothing to import any more
			if (n<1L){
				if (VERBOSE)
					cat("\n")
				break
			}
			## Recode levels if levels are specified
			if(any(appendLevels))
				for (i in i.fac){
					lev <- unique(c(levels(x[[i]]),levels(dat[[i]])))  # we save a call to the more general appendLevels() here
					levels(x[[i]]) <- lev
					dat[[i]] <- recodeLevels(dat[[i]], lev)
				}
			##
			## Append to the ffdf object
			##
			nff <- nrow(x)
			nrow(x) <- nff + n
			i <- hi(nff+1L, nff+n)
			x[i,] <- dat
			
			if (VERBOSE){
				write.stop <- proc.time()[3]
				write.time <- write.time + (write.stop - write.start)
				cat(" ffdf-write=", round(write.stop-write.start, 3), "sec\n", sep="")
			}			
			if(n<next.nrows){				
				break
			}				
		}		
	}
	if (VERBOSE){
		cat(" dbi-read=", round(read.time, 3), "sec  ffdf-write=", round(write.time, 3), "sec  TOTAL=", round(read.time+write.time, 3), "sec\n", sep="")
	}	
	return(x)
}


#' Read data from a ODBC connection into an ffdf.
#' 
#' Read data from a ODBC connection into an \code{\link[ff]{ffdf}}. This can for example be used to import
#' large datasets from Oracle, SQLite, MySQL, PostgreSQL, Hive or other SQL databases into R. \cr
#' 
#' Opens up the ODBC connection using \code{RODBC::odbcConnect}, sends the query using \code{RODBC::odbcQuery} and retrieves
#' the results in batches of next.rows rows using \code{RODBC::sqlGetResults}. Heavily borrowed from \code{\link[ff]{read.table.ffdf}}
#'
#' @param query the SQL query to execute on the ODBC connection
#' @param odbcConnect.args a list of arguments to pass to ODBC's \code{\link[RODBC]{odbcConnect}} (like dsn, uid, pwd). See the examples.
#' @param odbcQuery.args a list of arguments to pass to ODBC's \code{\link[RODBC]{odbcQuery}}, like rows_at_time. Defaults to an empty list.  
#' @param sqlGetResults.args a list containing optional parameters which will be passed to \code{\link[RODBC]{sqlGetResults}}.
#' Defaults to an empty list. The max parameter will be overwritten with first.rows and next.rows when importing in batches.     
#' @param x NULL or an optional ffdf object to which the read records are appended. 
#' See documentation in \code{\link[ff]{read.table.ffdf}} for more details and the example below.
#' @param nrows Number of rows to read from the query resultset. Default value of -1 reads in all rows.
#' @param first.rows chunk size (rows) to read for first chunk from the query resultset
#' @param next.rows chunk size (rows) to read sequentially for subsequent chunks from the query resultset. Currently, this must be specified.
#' @param levels optional specification of factor levels. A list with as names the names the columns of the data.frame 
#' fetched in the first.rows, containing levels of the factors. 
#' @param appendLevels logical. A vector of permissions to expand levels for factor columns. See documentation in \code{\link[ff]{read.table.ffdf}} for more details.
#' @param asffdf_args further arguments passed to \code{\link[ff]{as.ffdf}} (ignored if 'x' gives an ffdf object)
#' @param BATCHBYTES integer: bytes allowed for the size of the data.frame storing the result of reading one chunk. 
#' See documentation in \code{\link[ff]{read.table.ffdf}} for more details.
#' @param VERBOSE logical: TRUE to verbose timings for each processed chunk (default FALSE).
#' @param colClasses See documentation in \code{\link[ff]{read.table.ffdf}}
#' @param transFUN function applied to the data frame after each chunk is retreived by \code{\link[RODBC]{sqlGetResults}}
#' @param ... optional parameters passed on to transFUN
#' @return 
#' An ffdf object unless the query returns zero records in which case the function will return the data.frame
#' returned by \code{\link[RODBC]{sqlGetResults}} and possibly transFUN. 
#' @export
#' @seealso \code{\link[ff]{read.table.ffdf}, \link{read.dbi.ffdf}}
#' @examples
#' ##
#' ## Using the sqlite database (smalldb.sqlite3) in the /inst folder of the package
#' ## set up the sqlite ODBC driver (www.stats.ox.ac.uk/pub/bdr/RODBC-manual.pd) 
#' ## and call it 'smalltestsqlitedb' 
#' ##
#' \dontrun{
#' require(RODBC)
#' x <- read.odbc.ffdf(
#' query = "select * from testdata limit 10000",
#' odbcConnect.args = list(
#'  dsn="smalltestsqlitedb", uid = "", pwd = "", 
#'  believeNRows = FALSE, rows_at_time = 1), 
#' nrows = -1, 
#' first.rows = 100, next.rows = 1000, VERBOSE = TRUE)
#' }
read.odbc.ffdf <- function(
  query = NULL,
  odbcConnect.args = list(dsn=NULL, uid = "", pwd = ""), 
  odbcQuery.args = list(),
  sqlGetResults.args = list(), 		
  x = NULL, nrows = -1, 
  first.rows = NULL, next.rows = NULL, levels = NULL, appendLevels = TRUE, 
  asffdf_args = list(), BATCHBYTES = getOption("ffbatchbytes"), VERBOSE = FALSE, colClasses = NULL, 
  transFUN = NULL, ...){
  
  odbcinfo <- list()
  
  ##
  ## Connect to database
  ##
  odbcinfo$channel <- do.call('odbcConnect', odbcConnect.args)
  on.exit(try(RODBC::odbcClose(odbcinfo$channel), silent = TRUE))
  
  append <- !is.null(x)
  if(append && !inherits(x, "ffdf")){
    stop("only ffdf objects can be used for appending (and skipping the first.row chunk)")
  }
  if(VERBOSE){
    read.time <- 0
    write.time <- 0
  }
  
  nrows <- as.integer(nrows)
  N <- 0L	
  if(!append){		
    if(VERBOSE){
      cat("read.odbc.ffdf ", N+1L, "..", sep="")
      read.start <- proc.time()[3]
    }		
    if(is.null(colClasses)){
      colClasses <- NA
    }		
    ## Specify the size of the first number of rows to import
    if (is.null(first.rows)){
      if (is.null(next.rows)){
        first.rows <- 1e3L
      }else{
        first.rows <- as.integer(next.rows)
      }				
    }else{
      first.rows <- as.integer(first.rows)
    }		
    if (nrows>=0L && nrows<first.rows){
      first.rows <- nrows
    }		
    if (first.rows==1){
      stop("first.rows must not be 1")
    }
    ##
    ## Set up arguments for odbcQuery 
    ##
    odbcQuery.args$channel <- odbcinfo$channel
    odbcQuery.args$query <- query
    ##
    ## Launch the query
    ##
    odbcinfo$resultset <- do.call('odbcQuery', odbcQuery.args)
    
    ##
    ## Read in the first chunk
    ##
    sqlGetResults.args$channel <- odbcQuery.args$channel
    sqlGetResults.args$max <- first.rows			
    if(is.null(transFUN)){
      dat <- convertCharacterToFactor(do.call("sqlGetResults", sqlGetResults.args))
    }else{
      dat <- convertCharacterToFactor(transFUN(do.call("sqlGetResults", sqlGetResults.args), ...))
    }				
    n <- nrow(dat)
    N <- n
    ## Recode levels if levels are specified
    if(!is.null(levels)){
      cnam <- colnames(dat)
      lnam <- names(levels)
      if(!is.list(levels) || is.null(lnam) || any(is.na(match(lnam,cnam)))){
        stop("levels must be a list with names matching column names of the first data.frame read")
      }				
      for (i in lnam){
        dat[[i]] <- recodeLevels(dat[[i]], levels[[i]])
      }
    }		
    if (VERBOSE){
      write.start <- proc.time()[3]
      read.time <- read.time + (write.start - read.start)
      cat(N, " (", n, ")  odbc-read=", round(write.start-read.start, 3), "sec", sep="")
    }
    ##
    ## If there are zero rows in the query, return a data.frame instead of an ffdf object
    ##
    if(nrow(dat) == 0){
      if (VERBOSE){
        cat(" query returned 0 records\n", sep="")
        cat(" odbc-read=", round(read.time, 3), "sec  TOTAL=", round(read.time, 3), "sec\n", sep="")
      }				
      return(dat)
    }
    
    ##
    ## If there are non-zero rows in the query, convert to a filebased ffdf object
    ##
    x <- do.call("as.ffdf", c(list(dat), asffdf_args))		
    ##
    ## Fix ordered factors
    ##
    colClasses <- repnam(colClasses, colnames(x), default=NA)
    i.fix <- seq.int(length.out=ncol(dat))[!is.na(match(colClasses, "ordered"))]
    for (i in i.fix){
      virtual(x[[i]])$ramclass <- c("ordered","factor")
    }		
    if (VERBOSE){
      write.stop <- proc.time()[3]
      write.time <- write.time + (write.stop - write.start)
      cat(" ffdf-write=", round(write.stop-write.start, 3), "sec\n", sep="")
    }
  }
  ##
  ## Fetch the other rows
  ##
  if (append || N==first.rows){
    ##
    ## Define the number of rows to fetch in each run depending on the BATCHBYTES maximum
    ##
    k <- ncol(x)
    col.names <- colnames(x)
    colClasses <- sapply(seq.int(length.out=ncol(x)), function(i)colClass(x[[i]]))		
    if (is.null(next.rows)){
      recordsize <- sum(.rambytes[vmode(x)])
      next.rows <- BATCHBYTES %/% recordsize
      if (next.rows<1L){
        next.rows <- 1L
        warning("single record does not fit into BATCHBYTES")
      }
    }else{
      next.rows <- as.integer(next.rows)
    }		
    next.nrows <- next.rows
    
    appendLevels <- repnam(appendLevels, col.names, default=TRUE)
    if(any(appendLevels)){
      i.fac <- seq.int(length.out=k)
      if(append == TRUE){
        i.fac <- i.fac[appendLevels & sapply(i.fac, function(i) is.factor(x[[i]]))]
      }else{
        i.fac <- i.fac[appendLevels & sapply(i.fac, function(i) is.factor(dat[[i]]))]	
      }
      
    }	
    if(append == TRUE){
      ## Set up arguments for odbcQuery 
      odbcQuery.args$channel <- odbcinfo$channel
      odbcQuery.args$query <- query
      odbcinfo$resultset <- do.call('odbcQuery', odbcQuery.args)
    }
    ##
    ## Loop - fetching the subsequent data in batches
    ##
    while(TRUE){
      if (nrows>=0L && N+next.rows > nrows){
        next.nrows <- nrows - N
      }				
      if (next.nrows<1L){
        break
      }			
      if (VERBOSE){
        cat("read.odbc.ffdf ", N+1L,"..", sep="")
        read.start <- proc.time()[3]
      }
      ##
      ## read in subsequent chunk
      ##
      sqlGetResults.args$channel <- odbcQuery.args$channel
      sqlGetResults.args$max <- next.nrows
      if (is.null(transFUN)){
        dat <- convertCharacterToFactor(do.call("sqlGetResults", sqlGetResults.args))
      }else{
        dat <- convertCharacterToFactor(transFUN(do.call("sqlGetResults", sqlGetResults.args), ...))
      }			
      n <- nrow(dat)
      N <- N + n
      if (VERBOSE){
        write.start <- proc.time()[3]
        read.time <- read.time + (write.start - read.start)
        cat(N, " (", n, ")  odbc-read=", round(write.start-read.start, 3), "sec", sep="")
      }
      ## Nothing to import any more
      if (n<1L){
        if (VERBOSE)
          cat("\n")
        break
      }
      ## Recode levels if levels are specified
      if(any(appendLevels))
        for (i in i.fac){
          lev <- unique(c(levels(x[[i]]),levels(dat[[i]])))  # we save a call to the more general appendLevels() here
          levels(x[[i]]) <- lev
          dat[[i]] <- recodeLevels(dat[[i]], lev)
        }
      ##
      ## Append to the ffdf object
      ##
      nff <- nrow(x)
      nrow(x) <- nff + n
      i <- hi(nff+1L, nff+n)
      x[i,] <- dat
      
      if (VERBOSE){
        write.stop <- proc.time()[3]
        write.time <- write.time + (write.stop - write.start)
        cat(" ffdf-write=", round(write.stop-write.start, 3), "sec\n", sep="")
      }			
      if(n<next.nrows){				
        break
      }				
    }		
  }
  if (VERBOSE){
    cat(" odbc-read=", round(read.time, 3), "sec  ffdf-write=", round(write.time, 3), "sec  TOTAL=", round(read.time+write.time, 3), "sec\n", sep="")
  }	
  return(x)
}

#' Read data from a JDBC connection into an ffdf.
#' 
#' Read data from a JDBC connection into an \code{\link[ff]{ffdf}}. This can for example be used to import
#' large datasets from Oracle, SQLite, MySQL, PostgreSQL, Hive or other SQL databases into R. \cr
#' 
#' Opens up the JDBC connection using \code{RJDBC::dbConnect}, sends the query using \code{RJDBC::dbSendQuery} and \code{RJDBC::fetch}-es 
#' the results in batches of next.rows rows. Heavily borrowed from \code{\link[ff]{read.table.ffdf}} 
#'
#' @param query the SQL query to execute on the JDBC connection
#' @param dbConnect.args a list of arguments to pass to JDBC's \code{RJDBC::dbConnect} (like drv, dbname, username, password). See the examples.
#' @param dbSendQuery.args a list containing database-specific parameters which will be passed to to pass to \code{RJDBC::dbSendQuery}.
#' Defaults to an empty list.  
#' @param fetch.args a list containing optional database-specific parameters which will be passed to to pass to \code{RJDBC::fetch}. 
#' Defaults to an empty list.          
#' @param x NULL or an optional ffdf object to which the read records are appended. 
#' See documentation in read.table.ffdf for more details and the example below.
#' @param nrows Number of rows to read from the query resultset. Default value of -1 reads in all rows.
#' @param first.rows chunk size (rows) to read for first chunk from the query resultset
#' @param next.rows chunk size (rows) to read sequentially for subsequent chunks from the query resultset. Currently, this must be specified.
#' @param levels optional specification of factor levels. A list with as names the names the columns of the data.frame 
#' fetched in the first.rows, containing levels of the factors. 
#' @param appendLevels logical. A vector of permissions to expand levels for factor columns. See documentation in \code{\link[ff]{read.table.ffdf}} for more details.
#' @param asffdf_args further arguments passed to \code{\link[ff]{as.ffdf}} (ignored if 'x' gives an ffdf object )
#' @param BATCHBYTES integer: bytes allowed for the size of the data.frame storing the result of reading one chunk. 
#' See documentation in \code{\link[ff]{read.table.ffdf}} for more details.
#' @param VERBOSE logical: TRUE to verbose timings for each processed chunk (default FALSE).
#' @param colClasses See documentation in \code{\link[ff]{read.table.ffdf}}
#' @param transFUN function applied to the data frame after each chunk is retreived by \code{RJDBC::fetch}
#' @param ... optional parameters passed on to transFUN
#' @return 
#' An ffdf object unless the query returns zero records in which case the function will return the data.frame
#' returned by \code{RJDBC::fetch} and possibly transFUN. 
#' @export
#' @seealso \code{\link[ff]{read.table.ffdf}, \link{read.jdbc.ffdf}}
#' @examples
#' \dontrun{
#' require(ff)
#' 
#' ##
#' ## Example query using data in sqlite
#' ##
#' require(RSQLite)
#' dbfile <- system.file("smalldb.sqlite3", package="ETLUtils")
#' drv <- JDBC(driverClass = "org.sqlite.JDBC", classPath = "/usr/local/lib/sqlite-jdbc-3.7.2.jar")
#' query <- "select * from testdata limit 10000"
#' x <- read.jdbc.ffdf(query = query, 
#'  dbConnect.args = list(drv = drv, url = sprintf("jdbc:sqlite:%s", dbfile)), 
#'  first.rows = 100, next.rows = 1000, VERBOSE=TRUE)
#' class(x)
#' x[1:10, ]
#' }
read.jdbc.ffdf <- function(
  query = NULL,
  dbConnect.args = list(drv=NULL, dbname = NULL, username = "", password = ""), 
  dbSendQuery.args = list(), 
  fetch.args = list(),   	
  x = NULL, nrows = -1, 
  first.rows = NULL, next.rows = NULL, levels = NULL, appendLevels = TRUE, 
  asffdf_args = list(), BATCHBYTES = getOption("ffbatchbytes"), VERBOSE = FALSE, colClasses = NULL, 
  transFUN = NULL, ...){

  cleanupConnection <- function(x){
    if("resultset" %in% names(x)){
      RJDBC::dbClearResult(res=x$resultset)
    }
    if("channel" %in% names(x)){
      RJDBC::dbDisconnect(x$channel)
    }
  }
  jdbcinfo <- list()
  
  ##
  ## Connect to database
  ##
  jdbcinfo$channel <- do.call('dbConnect', dbConnect.args)
  on.exit(cleanupConnection(jdbcinfo))
  
  append <- !is.null(x)
  if(append && !inherits(x, "ffdf")){
    stop("only ffdf objects can be used for appending (and skipping the first.row chunk)")
  }
  if(VERBOSE){
    read.time <- 0
    write.time <- 0
  }
  
  nrows <- as.integer(nrows)
  N <- 0L	
  if(!append){		
    if(VERBOSE){
      cat("read.jdbc.ffdf ", N+1L, "..", sep="")
      read.start <- proc.time()[3]
    }		
    if(is.null(colClasses)){
      colClasses <- NA
    }		
    ## Specify the size of the first number of rows to import
    if (is.null(first.rows)){
      if (is.null(next.rows)){
        first.rows <- 1e3L
      }else{
        first.rows <- as.integer(next.rows)
      }				
    }else{
      first.rows <- as.integer(first.rows)
    }		
    if (nrows>=0L && nrows<first.rows){
      first.rows <- nrows
    }		
    if (first.rows==1){
      stop("first.rows must not be 1")
    }
    ##
    ## Set up arguments for odbcQuery and sqlGetResults functions
    ##
    dbSendQuery.args$conn <- jdbcinfo$channel
    dbSendQuery.args$statement <- query
    ##
    ## Launch the query
    ##
    jdbcinfo$resultset <- do.call('dbSendQuery', dbSendQuery.args)
    
    ##
    ## Read in the first chunk
    ##
    fetch.args$res <- jdbcinfo$resultset
    fetch.args$n <- first.rows			
    if(is.null(transFUN)){
      dat <- convertCharacterToFactor(do.call("fetch", fetch.args))
    }else{
      dat <- convertCharacterToFactor(transFUN(do.call("fetch", fetch.args), ...))
    }				
    n <- nrow(dat)
    N <- n
    ## Recode levels if levels are specified
    if(!is.null(levels)){
      cnam <- colnames(dat)
      lnam <- names(levels)
      if(!is.list(levels) || is.null(lnam) || any(is.na(match(lnam,cnam)))){
        stop("levels must be a list with names matching column names of the first data.frame read")
      }				
      for (i in lnam){
        dat[[i]] <- recodeLevels(dat[[i]], levels[[i]])
      }
    }		
    if (VERBOSE){
      write.start <- proc.time()[3]
      read.time <- read.time + (write.start - read.start)
      cat(N, " (", n, ")  jdbc-read=", round(write.start-read.start, 3), "sec", sep="")
    }
    ##
    ## If there are zero rows in the query, return a data.frame instead of an ffdf object
    ##
    if(nrow(dat) == 0){
      if (VERBOSE){
        cat(" query returned 0 records\n", sep="")
        cat(" jdbc-read=", round(read.time, 3), "sec  TOTAL=", round(read.time, 3), "sec\n", sep="")
      }				
      return(dat)
    }
    
    ##
    ## If there are non-zero rows in the query, convert to a filebased ffdf object
    ##
    x <- do.call("as.ffdf", c(list(dat), asffdf_args))		
    ##
    ## Fix ordered factors
    ##
    colClasses <- repnam(colClasses, colnames(x), default=NA)
    i.fix <- seq.int(length.out=ncol(dat))[!is.na(match(colClasses, "ordered"))]
    for (i in i.fix){
      virtual(x[[i]])$ramclass <- c("ordered","factor")
    }		
    if (VERBOSE){
      write.stop <- proc.time()[3]
      write.time <- write.time + (write.stop - write.start)
      cat(" ffdf-write=", round(write.stop-write.start, 3), "sec\n", sep="")
    }
  }
  ##
  ## Fetch the other rows
  ##
  if (append || N==first.rows){
    ##
    ## Define the number of rows to fetch in each run depending on the BATCHBYTES maximum
    ##
    k <- ncol(x)
    col.names <- colnames(x)
    colClasses <- sapply(seq.int(length.out=ncol(x)), function(i)colClass(x[[i]]))		
    if (is.null(next.rows)){
      recordsize <- sum(.rambytes[vmode(x)])
      next.rows <- BATCHBYTES %/% recordsize
      if (next.rows<1L){
        next.rows <- 1L
        warning("single record does not fit into BATCHBYTES")
      }
    }else{
      next.rows <- as.integer(next.rows)
    }		
    next.nrows <- next.rows
    
    appendLevels <- repnam(appendLevels, col.names, default=TRUE)
    if(any(appendLevels)){
      i.fac <- seq.int(length.out=k)
      if(append == TRUE){
        i.fac <- i.fac[appendLevels & sapply(i.fac, function(i) is.factor(x[[i]]))]
      }else{
        i.fac <- i.fac[appendLevels & sapply(i.fac, function(i) is.factor(dat[[i]]))]	
      }
      
    }	
    if(append == TRUE){
      ## Set up arguments for odbcQuery and sqlGetResults functions & Launch the query
      dbSendQuery.args$conn <- jdbcinfo$channel
      dbSendQuery.args$statement <- query			
      jdbcinfo$resultset <- do.call('dbSendQuery', dbSendQuery.args)
    }
    ##
    ## Loop - fetching the subsequent data in batches
    ##
    while(TRUE){
      if (nrows>=0L && N+next.rows > nrows){
        next.nrows <- nrows - N
      }				
      if (next.nrows<1L){
        break
      }			
      if (VERBOSE){
        cat("read.jdbc.ffdf ", N+1L,"..", sep="")
        read.start <- proc.time()[3]
      }
      ##
      ## read in subsequent chunk
      ##
      fetch.args$res <- jdbcinfo$resultset
      fetch.args$n <- next.nrows
      if (is.null(transFUN)){
        dat <- convertCharacterToFactor(do.call("fetch", fetch.args))
      }else{
        dat <- convertCharacterToFactor(transFUN(do.call("fetch", fetch.args), ...))
      }			
      n <- nrow(dat)
      N <- N + n
      if (VERBOSE){
        write.start <- proc.time()[3]
        read.time <- read.time + (write.start - read.start)
        cat(N, " (", n, ")  jdbc-read=", round(write.start-read.start, 3), "sec", sep="")
      }
      ## Nothing to import any more
      if (n<1L){
        if (VERBOSE)
          cat("\n")
        break
      }
      ## Recode levels if levels are specified
      if(any(appendLevels))
        for (i in i.fac){
          lev <- unique(c(levels(x[[i]]),levels(dat[[i]])))  # we save a call to the more general appendLevels() here
          levels(x[[i]]) <- lev
          dat[[i]] <- recodeLevels(dat[[i]], lev)
        }
      ##
      ## Append to the ffdf object
      ##
      nff <- nrow(x)
      nrow(x) <- nff + n
      i <- hi(nff+1L, nff+n)
      x[i,] <- dat
      
      if (VERBOSE){
        write.stop <- proc.time()[3]
        write.time <- write.time + (write.stop - write.start)
        cat(" ffdf-write=", round(write.stop-write.start, 3), "sec\n", sep="")
      }			
      if(n<next.nrows){				
        break
      }				
    }		
  }
  if (VERBOSE){
    cat(" jdbc-read=", round(read.time, 3), "sec  ffdf-write=", round(write.time, 3), "sec  TOTAL=", round(read.time+write.time, 3), "sec\n", sep="")
  }	
  return(x)
}





#' Write ffdf data to a database table by using a DBI connection.
#' 
#' Write \code{\link[ff]{ffdf}} data to a database table by using a DBI connection.
#' This can for example be used to store large ffdf datasets from R in
#' Oracle, SQLite, MySQL, PostgreSQL, Hive or other SQL databases. \cr
#' Mark that for very large datasets, these SQL databases might have tools to speed up by bulk loading.
#' You might also consider that as an alternative to using this procedure.
#' 
#' Opens up the DBI connection using \code{DBI::dbConnect}, writes data to the SQL table
#' using \code{DBI::dbWriteTable} by extracting the data in batches from the \code{\link[ff]{ffdf}}
#' and appending them to the table.
#'
#' @param x the \code{\link[ff]{ffdf}} to write to the database
#' @param name character string with the name of the table to store the data in. Passed on to \code{\link[DBI]{dbWriteTable}}.
#' @param dbConnect.args a list of arguments to pass to DBI's \code{\link[DBI]{dbConnect}} (like drv, dbname, username, password). See the examples.
#' @param BATCHBYTES integer: bytes allowed for the size of the data.frame storing the result of reading one chunk. 
#' See documentation in \code{\link[ff]{read.table.ffdf}} for more details.
#' @param RECORDBYTES optional integer scalar representing the bytes needed to process a single row of the ffdf
#' @param by integer passed on to \code{\link[bit]{chunk}} indicating to write to the database in chunks of this size. Overwrites
#' the behaviour of BATCHBYTES and RECORDBYTES.
#' @param VERBOSE logical: TRUE to verbose timings for each processed chunk (default FALSE).
#' @param ... optional parameters passed on to \code{\link[DBI]{dbWriteTable}}
#' @return 
#' invisible()
#' @export
#' @seealso \code{\link[DBI]{dbWriteTable}}, \code{\link[bit]{chunk}}
#' @examples
#' require(ff)
#' 
#' ##
#' ## Example query using data in sqlite
#' ##
#' require(RSQLite)
#' dbfile <- system.file("smalldb.sqlite3", package="ETLUtils")
#' drv <- dbDriver("SQLite")
#' query <- "select * from testdata limit 10000"
#' x <- read.dbi.ffdf(query = query, dbConnect.args = list(drv = drv, dbname = dbfile), 
#' first.rows = 100, next.rows = 1000, VERBOSE=TRUE)
#' 
#' write.dbi.ffdf(x = x, name = "helloworld", row.names = FALSE, overwrite = TRUE,
#'   dbConnect.args = list(drv = drv, dbname = dbfile), 
#'   by = 1000, VERBOSE=TRUE)
#' test <- read.dbi.ffdf(query = "select * from helloworld", 
#'   dbConnect.args = list(drv = drv, dbname = dbfile))
#' 
#' \dontrun{
#' require(ROracle)
#' write.dbi.ffdf(x = x, name = "hellooracle", row.names = FALSE, overwrite = TRUE,
#'   dbConnect.args = list(drv = dbDriver("Oracle"), 
#'                         user = "YourUser", password = "YourPassword", dbname = "Mydatabase"), 
#'   VERBOSE=TRUE)
#' }
write.dbi.ffdf <- function(x, name, 
                           dbConnect.args = list(drv=NULL, dbname = NULL, username = "", password = ""),
                           RECORDBYTES = sum(.rambytes[vmode(x)]), 
                           BATCHBYTES = getOption("ffbatchbytes"),
                           by = NULL,
                           VERBOSE = FALSE,
                           ...){
  stopifnot(inherits(x, "ffdf"))
  stopifnot(nrow(x) > 0)

  cleanupConnection <- function(x){
    if("channel" %in% names(x)){
      DBI::dbDisconnect(x$channel)
    }
  }
  dbiinfo <- list()
  dbiinfo$channel <- do.call('dbConnect', dbConnect.args)
  on.exit(cleanupConnection(dbiinfo))
    
  chunks <- ff::chunk.ffdf(x, RECORDBYTES=RECORDBYTES, BATCHBYTES=BATCHBYTES, by=by)
  for(i in seq_along(chunks)){
    if (VERBOSE){
      cat(sprintf("%s dbWriteTable chunk %s/%s\n", Sys.time(), i, length(chunks)))
    } 
    chunkidx <- chunks[[i]]
    dbWriteTable.args <- list(...)
    dbWriteTable.args$conn <- dbiinfo$channel
    dbWriteTable.args$name <- name
    dbWriteTable.args$value <- x[chunkidx, ]
    if(i > 1 && "overwrite" %in% names(dbWriteTable.args)){
      dbWriteTable.args$overwrite <- FALSE
    }
    if(i > 1){
      dbWriteTable.args$append <- TRUE
    }
    do.call('dbWriteTable', dbWriteTable.args)
  }
  invisible()
}


#' Write ffdf data to a database table by using a JDBC connection.
#' 
#' Write \code{\link[ff]{ffdf}} data to a database table by using a JDBC connection.
#' This can for example be used to store large ffdf datasets from R in
#' Oracle, SQLite, MySQL, PostgreSQL, Hive or other SQL databases. \cr
#' Mark that for very large datasets, these SQL databases might have tools to speed up by bulk loading.
#' You might also consider that as an alternative to using this procedure.
#' 
#' Opens up the JDBC connection using \code{RJDBC::dbConnect}, writes data to the SQL table
#' using \code{RJDBC::dbWriteTable} by extracting the data in batches from the \code{\link[ff]{ffdf}}
#' and appending them to the table.
#'
#' @param x the \code{\link[ff]{ffdf}} to write to the database
#' @param name character string with the name of the table to store the data in. Passed on to \code{dbWriteTable}.
#' @param dbConnect.args a list of arguments to pass to JDBC's \code{RJDBC::dbConnect} (like drv, dbname, username, password). See the examples.
#' @param BATCHBYTES integer: bytes allowed for the size of the data.frame storing the result of reading one chunk. 
#' See documentation in \code{\link[ff]{read.table.ffdf}} for more details.
#' @param RECORDBYTES optional integer scalar representing the bytes needed to process a single row of the ffdf
#' @param by integer passed on to \code{\link[bit]{chunk}} indicating to write to the database in chunks of this size. Overwrites
#' the behaviour of BATCHBYTES and RECORDBYTES.
#' @param VERBOSE logical: TRUE to verbose timings for each processed chunk (default FALSE).
#' @param ... optional parameters passed on to \code{dbWriteTable}
#' @return 
#' invisible()
#' @export
#' @seealso \code{\link[RJDBC]{JDBCConnection-methods}}, \code{\link[bit]{chunk}}
#' @examples
#' \dontrun{
#' require(ff)
#' 
#' ##
#' ## Example query using data in sqlite
#' ##
#' require(RJDBC)
#' dbfile <- system.file("smalldb.sqlite3", package="ETLUtils")
#' drv <- JDBC(driverClass = "org.sqlite.JDBC", classPath = "/usr/local/lib/sqlite-jdbc-3.7.2.jar")
#' query <- "select * from testdata limit 10000"
#' x <- read.jdbc.ffdf(query = query, 
#'  dbConnect.args = list(drv = drv, url = sprintf("jdbc:sqlite:%s", dbfile)), 
#'  first.rows = 100, next.rows = 1000, VERBOSE=TRUE)
#'  
#' write.jdbc.ffdf(x = x, name = "helloworld", row.names = FALSE, overwrite = TRUE,
#'   dbConnect.args = list(drv = drv, url = sprintf("jdbc:sqlite:%s", dbfile)), 
#'   by = 1000, VERBOSE=TRUE)
#' test <- read.jdbc.ffdf(query = "select * from helloworld", 
#'   dbConnect.args = list(drv = drv, url = sprintf("jdbc:sqlite:%s", dbfile)))
#' }
write.jdbc.ffdf <- function(x, name, 
                            dbConnect.args = list(drv=NULL, dbname = NULL, username = "", password = ""),
                            RECORDBYTES = sum(.rambytes[vmode(x)]), 
                            BATCHBYTES = getOption("ffbatchbytes"),
                            by = NULL,
                            VERBOSE = FALSE,
                            ...){
  stopifnot(inherits(x, "ffdf"))
  stopifnot(nrow(x) > 0)

  cleanupConnection <- function(x){
    if("channel" %in% names(x)){
      RJDBC::dbDisconnect(x$channel)
    }
  }
  dbiinfo <- list()
  dbiinfo$channel <- do.call('dbConnect', dbConnect.args)
  on.exit(cleanupConnection(dbiinfo))
  
  chunks <- ff::chunk.ffdf(x, RECORDBYTES=RECORDBYTES, BATCHBYTES=BATCHBYTES, by=by)
  for(i in seq_along(chunks)){
    if (VERBOSE){
      cat(sprintf("%s dbWriteTable chunk %s/%s\n", Sys.time(), i, length(chunks)))
    } 
    chunkidx <- chunks[[i]]
    dbWriteTable.args <- list(...)
    dbWriteTable.args$conn <- dbiinfo$channel
    dbWriteTable.args$name <- name
    dbWriteTable.args$value <- x[chunkidx, ]
    if(i > 1 && "overwrite" %in% names(dbWriteTable.args)){
      dbWriteTable.args$overwrite <- FALSE
    }
    if(i > 1){
      dbWriteTable.args$append <- TRUE
    }
    do.call('dbWriteTable', dbWriteTable.args)
  }
  invisible()
}



#' Write ffdf data to a database table by using a ODBC connection.
#' 
#' Write \code{\link[ff]{ffdf}} data to a database table by using a ODBC connection.
#' This can for example be used to store large ffdf datasets from R in
#' Oracle, SQLite, MySQL, PostgreSQL, Hive or other SQL databases. \cr
#' Mark that for very large datasets, these SQL databases might have tools to speed up by bulk loading.
#' You might also consider that as an alternative to using this procedure.
#' 
#' Opens up the ODBC connection using \code{RODBC::odbcConnect}, writes data to the SQL table
#' using \code{RODBC::sqlSave} by extracting the data in batches from the \code{\link[ff]{ffdf}}
#' and appending them to the table.
#'
#' @param x the \code{\link[ff]{ffdf}} to write to the database
#' @param tablename character string with the name of the table to store the data in. Passed on to \code{\link[RODBC]{sqlSave}}.
#' @param odbcConnect.args a list of arguments to pass to ODBC's \code{\link[RODBC]{odbcConnect}} (like dsn, uid, pwd). See the examples.
#' @param BATCHBYTES integer: bytes allowed for the size of the data.frame storing the result of reading one chunk. 
#' See documentation in \code{\link[ff]{read.table.ffdf}} for more details.
#' @param RECORDBYTES optional integer scalar representing the bytes needed to process a single row of the ffdf
#' @param by integer passed on to \code{\link[bit]{chunk}} indicating to write to the database in chunks of this size. Overwrites
#' the behaviour of BATCHBYTES and RECORDBYTES.
#' @param VERBOSE logical: TRUE to verbose timings for each processed chunk (default FALSE).
#' @param ... optional parameters passed on to \code{\link[RODBC]{sqlSave}}
#' @return 
#' invisible()
#' @export
#' @seealso \code{\link[RODBC]{sqlSave}}, \code{\link[bit]{chunk}}
#' @examples
#' ##
#' ## Using the sqlite database (smalldb.sqlite3) in the /inst folder of the package
#' ## set up the sqlite ODBC driver (www.stats.ox.ac.uk/pub/bdr/RODBC-manual.pd) 
#' ## and call it 'smalltestsqlitedb' 
#' ##
#' \dontrun{
#' require(RODBC)
#' x <- read.odbc.ffdf(
#'   query = "select * from testdata limit 10000",
#'   odbcConnect.args = list(
#'    dsn="smalltestsqlitedb", uid = "", pwd = "", 
#'    believeNRows = FALSE, rows_at_time = 1), 
#'   nrows = -1, 
#'   first.rows = 100, next.rows = 1000, VERBOSE = TRUE)
#'   
#' write.odbc.ffdf(x = x, tablename = "testdata", rownames = FALSE, append = TRUE,
#'   odbcConnect.args = list(
#'    dsn="smalltestsqlitedb", uid = "", pwd = "", 
#'    believeNRows = FALSE, rows_at_time = 1),  
#'   by = 1000, VERBOSE=TRUE)
#' }
write.odbc.ffdf <- function(x, tablename, 
                            odbcConnect.args = list(dsn=NULL, uid = "", pwd = ""), 
                            RECORDBYTES = sum(.rambytes[vmode(x)]), 
                            BATCHBYTES = getOption("ffbatchbytes"),
                            by = NULL,
                            VERBOSE = FALSE,
                            ...){
  stopifnot(inherits(x, "ffdf"))
  stopifnot(nrow(x) > 0)

  odbcinfo <- list()
  odbcinfo$channel <- do.call('odbcConnect', odbcConnect.args)
  on.exit(try(RODBC::odbcClose(odbcinfo$channel), silent = TRUE))
  
  chunks <- ff::chunk.ffdf(x, RECORDBYTES=RECORDBYTES, BATCHBYTES=BATCHBYTES, by=by)
  for(i in seq_along(chunks)){
    if (VERBOSE){
      cat(sprintf("%s sqlSave chunk %s/%s\n", Sys.time(), i, length(chunks)))
    } 
    chunkidx <- chunks[[i]]
    sqlSave.args <- list(...)
    sqlSave.args$channel <- odbcinfo$channel
    sqlSave.args$tablename <- tablename
    sqlSave.args$dat <- x[chunkidx, ]
    if(i > 1){
      sqlSave.args$append <- TRUE
    }
    do.call('sqlSave', sqlSave.args)
  }
  invisible()
}


# NOTE that this colClass implementation works only because accidentally the last position of the oldClasses is needed
# for c("ordered","factor") read.table does not want "ordered"
# for c("POSIXt","POSIXct") "POSIXct" is in the last position only due to "historical error" (Gabor Grothendieck, r-help, 26.9.2009)  
colClass <- function(x){
  UseMethod("colClass")
}	
colClass.default <- function(x){
  cl <- class(x)
  cl[length(cl)]
}	
colClass.ff <- function(x){
  if (length(x))
    x <- x[1]
  else
    x <- x[]
  NextMethod()
}
convertCharacterToFactor <- function(x){
  chartofactor <- sapply(x, class)
  chartofactor <- chartofactor[chartofactor == "character"]
  for(convertme in names(chartofactor)){
    x[[convertme]] <- factor(x[[convertme]]) 
  }
  x
}
