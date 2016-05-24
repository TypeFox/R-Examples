#' Set path in a file name
#'
#' This function sets the path in a file name.
#' 
#' @param x     a \code{character} or a \code{filename}
#' @param path  new path to file
#' @return modified object of the original type
#' @export
#' 
#' @examples
#' x <- "path/data_norm.txt"
#' print(set_fpath(x, "new_path"))
#'
set_fpath <- function(x, path) UseMethod("set_fpath");

#' @export
set_fpath.filename <- function(x, path) {
	if (length(path) == 1) {
		# path is not already split, so split it
		path <- strsplit(path, .Platform$file.sep, fixed=TRUE)[[1]];
	}
	x$path <- path;
	x
}

#' @export
set_fpath.character <- function(x, path) {
	as.character(set_fpath.filename(as.filename(x), path))
}

#' Set file extension
#'
#' This function sets the extension in a file name.
#' 
#' @param x     a \code{character} or a \code{filename}
#' @param ext   new file extension
#' @return modified object of the original type
#' @export
#' 
#' @examples
#' x <- "data_norm_2011-01-03.txt"
#' print(set_fext(x, "csv"))
#'
set_fext <- function(x, ext) UseMethod("set_fext");

#' @export
set_fext.filename <- function(x, ext) {
	if (length(ext) == 1) {
		# extension is not already split, so split it
		ext <- strsplit(ext, .get_ext_char(), fixed=TRUE)[[1]];
	}
	x$ext <- ext;
	x
}

#' @export
set_fext.character <- function(x, ext) {
	as.character(set_fext.filename(as.filename(x), ext))
}

#' Set date stamp in a file name
#'
#' This function sets the date stamp in a file name.
#' 
#' @param x     a \code{character} or a \code{filename}
#' @param date  new date stamp (\code{character} or \code{Date})
#' @return modified object of the original type
#' @export
#' 
#' @examples
#' x <- "data_norm_2011-01-03.txt"
#' print(set_fdate(x, "2011-01-05"))
#'
set_fdate <- function(x, date) UseMethod("set_fdate");

#' @export
set_fdate.filename <- function(x, date) {
	if (is.character(date)) {
		if (!.grepl_date(date)) {
			stop("Invalid date format");
		}
		x$date <- date;
	} else if (is(date, "Date")) {
		x$date <- .format_date(date);
	} else {
		stop("`date` must be a character or a Date object");
	}
	x
}

#' @export
set_fdate.character <- function(x, date) {
	as.character(set_fdate.filename(as.filename(x), date))
}

#' Set time stamp in a file name
#'
#' This function sets the time stamp in a file name.
#' 
#' @param x     a \code{character} or a \code{filename}
#' @param time  new time stamp (\code{character} or \code{POSIXct})
#' @return modified object of the original type
#' @export
#' 
#' @examples
#' x <- "data_norm_20110103T093015.txt"
#' # change the time to 30 seconds past 2:45 p.m. 
#' print(set_ftime(x, "144530"))
#' # to change the date, time must be specified as well
#' print(set_ftime(x, "20110505T101500"))
#'
set_ftime <- function(x, time) UseMethod("set_ftime");

#' @export
set_ftime.filename <- function(x, time) {
	if (is.character(time)) {
		if (.grepl_time(time)) {
			if (is.na(x$date)) {
				stop("filename `x` has no date stamp; cannot apply only time stamp without date stamp");
			} else {
				x$time <- time;
			}
		} else if (.grepl_datetime(time)) {
			x <- .set_fdatetime(x, time);
		} else {
			stop("Invalid time or date-time format");
		}
	} else if (is(time, "POSIXct")) {
		x$date <- .format_date(time, compact=TRUE);
		x$time <- .format_time(time);
	} else {
		stop("`time` must be a character or POSIXct object");
	}
	x
}

#' @export
set_ftime.character <- function(x, time) {
	as.character(set_ftime.filename(as.filename(x), time))
}

