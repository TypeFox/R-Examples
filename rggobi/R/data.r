# Set data
# ===================================================================

# [<-.GGobi
# Add data to ggobi instance.
# 
# This function allows you to add (and eventually) replace
# GGobiData objects in a GGobi instance.  
# 
# @arguments ggobi instance
# @arguments name of data frame
# @arguments data.frame, or string to path of file to load
# @alias $<-.GGobi
# @alias [[<-.GGobi
# @keyword manip 
#X g <- ggobi()
#X g["a"] <- mtcars
#X g$b <- mtcars
"[<-.GGobi" <- function(x, i, value) {
	if (inherits(value, "GGobiData")) return(x)
	if (length(i) != 1) stop("You may only add or replace one dataset at a time")
	replace <- (i %in% names(x))
	if (replace) {
	  warning("Replacement does not work")
	  return(x)
	}
	
	if(is.character(value)) {
		ggobi_set_data_file(value, .gobi = x)
	} else {
		ggobi_set_data_frame(as.data.frame(value), name = clean.ggobi(i), .gobi = x)
	}
	x
}
"[[<-.GGobi" <- "$<-.GGobi" <- function(x, i, value) {
  x[i] <- value
  x
}
# Set data file.
# Open data file and add to ggobi datasets.
# 
# @arguments path to file
# @arguments mode of file
# @arguments add?
# @arguments ggobi instance
# @returns GGobiData
# @keyword manip 
# @keyword internal
ggobi_set_data_file <- function(file, mode = "unknown", add = TRUE, .gobi = ggobi_get()) {
	num <- .GGobiCall("setFile", as.character(file), as.character(mode), as.logical(add), .gobi = .gobi)
	if(num > -1) {
		dataset(num, .gobi)
	} else {
		warning(paste("Failed to open", file))
	}
}


# Set data frame.
# Add data.frame to ggobi instance.
# 
# @arguments data frame to add
# @arguments data set name (appears on tabs in ggobi)
# @arguments description of data frame
# @arguments rownames
# @arguments ggobi instance
# @returns GGobiData
# @keyword manip 
# @keyword internal
ggobi_set_data_frame <- function(data, name = deparse(sys.call()[[2]]), description = paste("R data frame", name), id = NULL, .gobi = ggobi_get()) {
	n <- dimnames(data)
	
	if (is.null(id)) {
	  if (.row_names_info(data) < 0) {
	    id <- paste(name, seq(1, nrow(data)), sep="")
	  } else {
	    id <- rownames(data)
	  }
	} else {
	  if (length(id) != nrow(data)) stop("Length of id does not match rows of data")
	}
	
	# Convert character columns to factors
	char <- unlist(lapply(data, is.character))
	data[char] <- lapply(data[char], factor)

	.data <- .GGobiCall("addData",
		data,
		n[[1]], n[[2]], dim(data), as.character(description),
		as.character(name),
		as.character(id), .gobi = .gobi)

	.data
}


# Get data
# ===================================================================

# Get ggobi data.
# Conveniently retrieve ggobi dataset.
#
# It is convenient to be able to refer to and operate on a ggobi 
# dataset as if it were a regular R dataset.  This function allows one to
# get an \code{GGobiData} object that represents a particular 
# dataset. 
# 
# @arguments GGobi object
# @arguments name of dataset to retrive
# @arguments ignored
# @arguments if TRUE, return vector is possible
# @keyword manip 
# @alias [[.GGobi
# @alias $.GGobi
#X g <- ggobi(ChickWeight)
#X g["cars"] <- mtcars
#X g[1:2]
#X g["ChickWeight"]
#X g["cars"]
#X g$cars
"[.GGobi" <- function(x, i, ..., drop=TRUE) {
	d <- dataset(clean.ggobi(i), .gobi = x)
	if (drop && length(d) == 1 ) {
		d[[1]]
	} else {
		d
	}
}
"[[.GGobi" <- "$.GGobi" <- function(x, i) {
  x[i]
}

# Generic method for getting dataset
# @keyword internal 
dataset <- function(x, .gobi = ggobi_get()) UseMethod("dataset", x)

# Get ggobi dataset.
# Get an object representing an internal ggobi dataset
#
# It is convenient to be able to refer to and operate on a ggobi 
# dataset as if it were a regular R dataset.  This function allows one to
# get an \code{GGobiData} object that represents a particular 
# dataset. 
# 
# @arguments which dataset to retrieve, an integer for positional matching or a character to match by name
# @arguments GGobi instance
# @value A list of \code{GGobiData} objects
# @seealso \code{\link{$.GGobiData}} for user level selection of datasets
# @keyword manip 
# @keyword internal
# @alias dataset.character
dataset.numeric <- function(x, .gobi = ggobi_get()) {
	refs <- .GGobiCall("getDataset", as.integer(x-1), .gobi=.gobi)

	refs <- lapply(refs, function(x) {
		if(is.null(x)) return() 
		#class(x) <- c(class(x), "data.frame"); 
		attr(x,"ggobi") <- .gobi
		x
	}) 

	refs
}
dataset.character <- function(x, .gobi = ggobi_get()) {
	id <- match(x, names(.gobi))
	if (any(is.na(id))) {
		stop(paste("Unrecognized dataset name", x[is.na(id)]))
	}
	
	dataset(id, .gobi)
}

# Write xml
# Write GGobi xml for specific dataset to filename
#
# @arguments GGobiData object
# @arguments path to write file to
# @keyword manip
# ggobi_data_write_xml <- function (gd, filename) {
# 	refs <- lapply(.data, dataset, .gobi)
# 	.GGobiCall("writeDatasetsXML", gd, as.character(filename))
# }
