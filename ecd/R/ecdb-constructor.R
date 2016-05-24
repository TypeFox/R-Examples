#' Constructor of ecdb class for the elliptic database
#' 
#' Construct an ecdb class by providing the required parameters.
#' The default is to use the internal database location. 
#' But the internal db is limited in size.
#' The the elliptic database stores
#' the stdev, kurtosis, discriminant, j-invariant, and ellipticity.
#' for alpha and gamma between -100 and 100.
#' Step size is 1 for -100 to 100; 0.25 for -50 to 50;
#' 0.1 for -10 to 10; 0.025 between -6 and 1.
#' Speical lines with step size of 0.001 for j0 and j1728 between -10 and 10;
#' 0.01 for kmax and critical between 0 and 100.
#' For asym1X, step size is 10 from 100 to 1000.
#' For asym2X, step size is 100 from 1000 to 10000.
#' For asym3X, step size is 1000 from 10000 to 60000.
#' For polar-q1, step size is 0.025 from 0 to 20 for log2(R), and integer angles, 0-89.
#'
#' @param file Character, the full path to an elliptic database. 
#'              Use "internal" to force the usage of the internal db.
#' @param newdb Logical. If \code{TRUE}, remove existing db and create a new one.
#'              Default: \code{FALSE}.
#'
#' @return An object of ecdb class 
#'
#' @keywords constructor
#'
#' @export ecdb
#' 
#' @examples
#' db <- ecdb("internal")

### <======================================================================>
"ecdb" <- function(file = NULL, newdb=FALSE)
{
	is.internal <- FALSE

	get_internal_file <- function() {
	    is.internal <<- TRUE
	    paste(system.file("extdata", package = "ecd"), "elliptic.db", sep="/")
	}
    
    # locate the db file
    if (toString(file) == "internal") {
        file <- get_internal_file()
    }
    else if (is.null(file)) {
    	if (length(getOption("elliptic_db")) > 0) {
        	file <- getOption("elliptic_db")
		} else {
			file <- paste("~", "elliptic.db", sep="/")
			if (! file.exists(file)) {
			    file <- get_internal_file()
			}
		}
    }
    
    if (newdb) {
        if (file.exists(file)) {
        	file.remove(file)
        	warning(paste("The elliptic db file is removed and recreated:", file))
        }
    }
	else if (! file.exists(file)) {
	    stop(paste("Failed to locate existing elliptic db:", file))
	}
	    
    conn <- RSQLite::dbConnect(SQLite(), dbname=file)
    call <- match.call()
    db <- new("ecdb", call = call,
               file = unname(file),
               conn = conn,
               is.internal = is.internal,
               conf = .ecdb.conf())

    if (newdb) {        
        # generate the schema
        ecdb.dbSendQuery(db, 
            "CREATE TABLE ECDATTR (
             alpha_m INTEGER,
             gamma_m INTEGER,
             cusp INTEGER,
             stdev   REAL,
             kurtosis REAL,
             discr REAL,
             jinv REAL,
             ellipticity REAL,
             const REAL,
             time_stamp INTEGER,
             PRIMARY KEY (alpha_m, gamma_m)
            )")

        ecdb.dbSendQuery(db, 
            "CREATE TABLE LAMBDAATTR (
             lambda_m INTEGER,
             beta_m INTEGER,
             sigma_m INTEGER,
             stdev REAL,
             kurtosis REAL,
             const REAL,
             mnt1 REAL,
             mnt2 REAL,
             trunc_intg REAL,
             trunc_mnt REAL,
             mlog_mgf1 REAL,
             time_stamp INTEGER,
             PRIMARY KEY (lambda_m, beta_m, sigma_m)
            )")
        
		# keep track of versions and upgrades
        ecdb.dbSendQuery(db, 
            "CREATE TABLE VERSION (
             version INTEGER PRIMARY KEY,
             note TEXT,
             iso_dttm TEXT,
             time_stamp INTEGER
            )")

        # keep track of write history
        ecdb.dbSendQuery(db,
            "CREATE TABLE HISTORY (
                iso_dttm TEXT PRIMARY KEY,
                activity TEXT,
                time_stamp INTEGER
            )")

		sql <- "INSERT INTO VERSION VALUES (:version, :note, :iso_dttm, :time_stamp)"
        version <- data.frame(version=c(104),
                              note=c("ecdb version 1.04: adds lambda attributes"),
                              iso_dttm=c(as.character(Sys.time())),
                              time_stamp=c(as.integer(Sys.time()))
                             )
    	rs <- RSQLite::dbGetPreparedQuery(conn, sql, bind.data=version)
        ecdb.protectiveCommit(db)

        history(db) <- "newdb"

        if (! file.exists(file)) {
        	stop(paste("Failed to create elliptic db at", file))
		}
	}
    
    return(db)
}
### <---------------------------------------------------------------------->
# internal helper to render the standard configuration
.ecdb.conf <- function() {
    
    # retangle of data on (a,r) plane
    # type = c(from, to, step)
    rect <- list(
        "100" = list(from= -100, to= 100, step=1),
        "50"  = list(from= -50,  to= 50,  step=0.25),
        "10"  = list(from= -10,  to= 10,  step=0.1),
        "6"   = list(from= -1,   to= 6,   step=0.025),
        # internal: this is all we can do to fit a 1MB sqlite db
        "internal" = list(from= -10, to= 10, step=0.25))

    # retangle of data on (R,theta) plane
    polar <- list(
        "q1"  = list( # polar, 1st quadrant
            log2R = list(from= 0, to= 20, step=0.025, var="log2R"),  # R=2^N
            angle = list(from= 0, to= 89, step=1,     var="angle")
            )
        #"q2"  = list( # polar, 2nd quadrant
        #    log2R = list(from= 0, to= 20, step=0.025, var="log2R"),  # R=2^N
        #    angle = list(from= 90, to= 179, step=1,   var="angle")
        #),
        #"q3"  = list( # polar, 3rd quadrant
        #    log2R = list(from= 0, to= 20, step=0.025, var="log2R"),  # R=2^N
        #    angle = list(from= 180, to= 269, step=1,  var="angle"), mpfr=1
        #),
        #"q4"  = list( # polar, 4th quadrant
        #    log2R = list(from= 0, to= 20, step=0.025, var="log2R"),  # R=2^N
        #    angle = list(from= 270, to= 314, step=1,  var="angle"), mpfr=1
        #)
    )
    
    # a specialized line of data
    # type = c(from, to, step, [var], [cusp])
    # var: default is alpha if not mentioned.
    # cusp: default is 0 if not mentioned
    line <- list(
        "j0"    = list(from= -10,  to= 10,  step=0.001, var="gamma", val=0),
        "j1728" = list(from= -10,  to= 10,  step=0.001, var="alpha", val=0),
        "kmax"  = list(from= -100, to= 100, step=0.01,  var="alpha", val=2.910),
        "critical" = list(from= 0, to= 100, step=0.01,  var="critical", cusp=1),
        "critical1" = list(from= 100, to= 1000, step=1,  var="critical", cusp=1, mpfr=1),
        "critical2" = list(from= 1000, to= 10000, step=10,  var="critical", cusp=1, mpfr=1),
        "critical3" = list(from= 10000, to= 60000, step=100,  var="critical", cusp=1, mpfr=1),
        
        # for asymptotic study, exponentially increasing steps
        # "asymp1" = list(from= 100,   to= 1000,  step=1,   var="4x", mpfr=1),
        # "asymp2" = list(from= 1000,  to= 10000, step=10,  var="4x", mpfr=1),
        # "asymp3" = list(from= 10000, to= 60000, step=100, var="4x", mpfr=1))
        # need more granularity to debug, so each 4x is broken up into 4 directions
        "asymp11" = list(from= 100,   to= 1000,  step=1,   var="gamma", val=0), # degree = 0
        "asymp12" = list(from= 100,   to= 1000,  step=1,   var="alpha", val=0), # degree = 90
        "asymp13" = list(from= -1000, to= -100,  step=1,   var="gamma", val=0, mpfr=1), # degree = 180
        "asymp14" = list(from= -1000, to= -100,  step=1,   var="alpha", val=0, mpfr=1), # degree = 270

        "asymp21" = list(from= 1000,   to= 10000,  step=10,   var="gamma", val=0), # degree = 0
        "asymp22" = list(from= 1000,   to= 10000,  step=10,   var="alpha", val=0), # degree = 90
        "asymp23" = list(from= -10000, to= -1000,  step=10,   var="gamma", val=0, mpfr=1), # degree = 180
        "asymp24" = list(from= -10000, to= -1000,  step=10,   var="alpha", val=0, mpfr=1), # degree = 270

        "asymp31" = list(from= 10000,   to= 60000,  step=100,   var="gamma", val=0), # degree = 0
        "asymp32" = list(from= 10000,   to= 60000,  step=100,   var="alpha", val=0), # degree = 90
        "asymp33" = list(from= -60000, to= -10000,  step=100,   var="gamma", val=0, mpfr=1), # degree = 180
        "asymp34" = list(from= -60000, to= -10000,  step=100,   var="alpha", val=0, mpfr=1) # degree = 270
    )
    
    # final list
    list(
        rect = rect,
        polar = polar,
        line = line)
}