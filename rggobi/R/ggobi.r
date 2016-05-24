ggobi <- function(data=NULL, ...) UseMethod("ggobi", data)

# New ggobi
# Creates a new ggobi instance
# 
# This function creates a new instance of GGobi with or without new data.  Use 
# this function whenever you want to create a new GGobi indepdent of the 
# others---they will not share linked plots.  If you want to add
# another dataset to an existing ggobi, please see \code{\link{[<-.GGobi}}
# 
# There are currently three basic types of functions that you
# can use with rggobi:
# 
# \itemize{
#   \item Data getting and setting: see \code{\link{[.GGobi}}, and \code{\link{[.GGobiData}}
#   \item "Automatic" brushing: see \code{\link{glyph_colour}}, 
#     \code{\link{glyph_size}},  \code{\link{glyph_type}},
#     \code{\link{shadowed}},    \code{\link{excluded}}, and the associated
#     setter functions.
#   \item Edge modifcation: see \code{\link{edges}}, \code{\link{edges<-}},
#     \code{\link{ggobi_longitudinal}}
# }
# 
# You will generally spend most of your time working with 
# \code{ggobdata}s, you retrieve using \code{\link{$.GGobiData}}, 
# \code{\link{[.GGobiData}}, or \code{\link{[[.GGobiData}}.
# Most of the time these will operate like normal R datasets while
# pointing to the data in GGobi so that all changes are kept in sync.  
# If you need to force a ggobiDaataset to be an R \code{data.frame} use
# \code{\link{as.data.frame}}.
#
# @arguments the name of a file containing the data, or a data frame or matrix containing the values
# @arguments a character vector of command-line arguments
# @arguments data format GGobi should expect to read the data from, if reading from a file.
# @arguments the name to use in GGobi for the dataset, if one is specified
# @arguments ignored
# @value A ggobi object 
# @keyword dynamic 
# @alias rggobi
# @alias ggobi
#X ggobi(ggobi_find_file("data", "flea.csv"))
#X ggobi(ggobi_find_file("data", "flea.xml"))
#X ggobi(mtcars)
#X mtcarsg <- ggobi_get()$mtcars
#X glyph_colour(mtcarsg)
#X glyph_colour(mtcarsg) <- ifelse(mtcarsg$cyl < 4, 1, 2)
#X glyph_size(mtcarsg) <- mtcarsg$cyl
ggobi.default <- function(data, args=character(0), mode=character(0), name = deparse(sys.call()[[2]]), ...) {
	
	filename <- character(0)
	if(!missing(data) && is.character(data) && file.exists(data)) {
		filename <- path.expand(data)
	}

	args <- c(file.path(getwd(),"rggobi"), "--keepalive", as.character(args), as.character(mode), filename)
  
	ok <- .Call(.ggobi.symbol("init"), args, TRUE, PACKAGE = "rggobi")

	if(!is.null(ok) && !missing(data) && length(filename) == 0) {
	  name <- force(name)
	  ok[name] <- data
	}

	invisible(ok)
}


# GGobi names
# Get dataset names
# 
# @arguments ggobi instance
# @keyword dynamic 
#X g <- ggobi(mtcars)
#X names(g)
names.GGobi <- function(x) {
 .GGobiCall("getDatasetNames", .gobi=x)
}

# Clean ggobi
# Clean arguments for ggobi
# 
# Arguments for ggobi need to be in specific format.
# This function helps ensure that.
# 
# @arguments vector
# @keyword dynamic 
# @keyword internal 
clean.ggobi <- function(x) {
	if (is.numeric(x)) {
		as.integer(x)
	} else {
		as.character(x)
	}
}


# Print ggobi
# Prints summary of ggobi object by instance
# 
# @arguments ggobi object
# @seealso \code{\link{summary.GGobi}}
# @keyword dynamic 
# @keyword internal 
print.GGobi <- function(x, ...) {
	print(summary(x))
}

# GGobi summary 
# Get a description of the global state of the GGobi session.
# 
# @arguments ggobi object
# @arguments ignored
# @keyword dynamic 
#X g <- ggobi(mtcars)
#X summary(g)
summary.GGobi <- function(object, ...) {
  ans <- .GGobiCall("getDescription", .gobi = object)
	if (is.null(ans)) return("Nothing known about this GGobi instance")
	
	dimnames(ans$"Data dimensions") <- list(ans$Fileame, c("nrow", "ncol"))
  ans
}

# Close GGobi instance
# Terminates and discards a ggobi instance
# 
# This allows the caller to close a ggobi instance and discard the
# resources it uses. The function closes the display windows and
# variable panel window associated with this ggobi instance.
# It also resets the default ggobi instance to be the last
# one created.
#
# @arguments ggobi object to close 
# @arguments ignored and for compatability generic function.
# @keyword dynamic 
#X g1 <- ggobi(mtcars)
#X g2 <- ggobi(mtcars)
#X close(g2)
#X close(ggobi_get())
close.GGobi <- function(con, ...) {
  ok <- .GGobiCall("close", .gobi = con)
  invisible(ok)
}

# Get number of GGobis
# Retrieves the number of ggobi instances within this session
#
# One or more ggobi instances can be created within an R session so that one 
# can simultaneously look at different datasets or have different views
# of the same dataset.  This function returns the number of ggobis in existence. 
# 
# The different ggobi instances are maintained as C  level structures. 
# This function accesses a variable that stores how many are in existence 
# when the function is invoked.
# 
# @keyword dynamic 
#X ggobi_count()
ggobi_count <- function() {
 .GGobiC("getNumGGobiInstances", num = as.integer(-1), .gobi=NULL)$num
}


# Get GGobi
# Returns a ggobi reference
#
# This allows one to get a list of all the ggobi instances currently
# in existence in the R session.  Also, one can fetch particular instances. 
# 
# @arguments numeric vector indicating which ggobi instances to retrieve.  Use default if none specified
# @arguments drop if possible?
# @returns list of ggobi instances
# @keyword dynamic 
#X ggobi(mtcars)
#X ggobi(Nile)
#X ggobi_get(1)
#X ggobi_get(1:2)
ggobi_get <- function(id = ggobi_count(), drop=TRUE) {
	ggobis <- .GGobiCall("getGGobi", as.integer(id), .gobi=NULL)

	if (drop && length(ggobis) == 1) return(ggobis[[1]])
	ggobis
}

# Get version
# GGobi version information
# 
# Return an object that describes the version of the ggobi 
# library being used. This allows code to execute 
# conditionally on certain version numbers, etc.
# 
# @value date the release date of the ggobi library
# @value version a vector of 3 integers containing the major, minor and patch-level numbers.
# @value versionstring a string version of the major, minor and patch-level numbers,
# @keyword dynamic 
#X ggobi_version()
ggobi_version <- function() {
	info <- .GGobiCall("getVersionInfo", .gobi=NULL)

	names(info) <- c("date", "version", "version string")
	names(info$version) <- c("major", "minor", "patch")
	info
}

# Find GGobi file.
# Finds a file stored somewhere in the ggobi installation.
# 
# @arguments bits of the path to join together
# @keyword dynamic 
# @keyword internal 
#X ggobi_find_file("data","tips.xml")
ggobi_find_file <- function(..., check = F) {
	file <- .GGobiCall("ggobi_find_data_file", file.path(...))
	if(check && !file.exists(file)) stop("Cannot find file ", file)
	
	file
}