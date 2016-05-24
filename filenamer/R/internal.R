# Get the default delimiter characters
.get_tag_char <- function(tag.char=NULL) {
	if (is.null(tag.char)) tag.char <- getOption("filename.tag.char");
	if (is.null(tag.char)) '_' else tag.char
}

# Get the default extension character
.get_ext_char <- function(ext.char=NULL) {
	if (is.null(ext.char)) ext.char <- getOption("filename.ext.char");
	if (is.null(ext.char)) "." else ext.char
}

# Format Date object
.format_date <- function(x, compact=FALSE) {
	if (compact) {
		format(x, "%Y%m%d");
	} else {
		format(x, "%Y-%m-%d");
	}
}

# Format POSIXct object
.format_time <- function(x) {
	format(x, "%H%M%S");
}

# Regular expression match against date format
.grepl_date <- function(x) {
	if (is.null(x)) {
		FALSE
	} else {
		grepl("^\\d{4}-\\d{2}-\\d{2}(-[A-Z]+)?$", x)
	}
}

# Regular expression match against date format
.grepl_time <- function(x) {
	if (is.null(x)) {
		FALSE
	} else {
		grepl("^\\d{6}$", x)
	}
}

# Regular expression match against date format
# Date and time are separated by "T", as per ISO 8601 specification
.grepl_datetime <- function(x) {
	if (is.null(x)) {
		FALSE
	} else {
		grepl("^\\d{4}-?\\d{2}-?\\d{2}T\\d{6}$", x)
	}
}

# Set datetime for a filename
# @param x         filename
# @param datetime  datetime character vector (check omitted)
.set_fdatetime <- function(x, datetime) {
		# date and time are separated by a 'T'
		datetime <- strsplit(datetime, "T", fixed=TRUE)[[1]];
		x$date <- datetime[1];
		x$time <- datetime[2];
		x
}

# Stamp a filename or filepath with date and time
# @param y     filename or filepath
# @param date  date character vector
# @param time  time character vector
.append_datetime <- function(y, date, time) {
	if (!is.na(date)) {
		if (!is.na(time)) {
			y <- c(y, paste(date, time, sep="T"));
		} else {
			y <- c(y, date);
		}
	}
	y
}

# Create a date-time stamp
# @param stamp  stamp option
.stamp_datetime <- function(date, time, stamp) {
	if (is.null(date) && is.null(time)) {
		# use automated stamping
		if (is.null(stamp) || stamp <= 0) {
			# no stamping
			date <- NA;
			time <- NA;
		} else if (stamp == 1) {
			# only date stamping
			date <- format(Sys.Date(), "%Y-%m-%d");
			time <- NA;
		} else {
			# date-time stamping
			s <- Sys.time();
			date <- format(s, "%Y%m%d");
			time <- format(s, "%H%M%S");
		}
	} else if (is(time, "POSIXct")) {
		# set both date and time stamp based on `time` parameter
		date <- format(time, "%Y%m%d");
		time <- format(time, "%H%M%S");
	} else if (is(date, "Date")) {
		# set only date stamp
		date <- as.character(date);
		time <- NA;
	} else if (is.null(time)) {
		# time alone is NULL, date is specified: no time stamping
		time <- NA;
	}

	# date and time should already have been converted to character or NA
	stopifnot(is.character(date) || is.na(date));
	stopifnot(is.character(time) || is.na(time));

	list(date=date, time=time)
}
