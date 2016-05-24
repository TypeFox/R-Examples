
fwf2csv <- function(fwffile, csvfile, names, begin, end, verbose = getOption("verbose"))
{
    # Check for errors
    ncols = length(names)
    if(length(begin) != ncols || length(end) != ncols)
	stop("The vectors \"names\", \"begin\" and \"end\" must have the same length.")
    if(file.exists(fwffile) == FALSE){
	msg <- paste(gettext("File not found:", domain = "R-descr"), fwffile)
	stop(msg)
    }

    csvfile <- path.expand(csvfile)
    fwffile <- path.expand(fwffile)

    .C("realfwf2csv",
	as.character(fwffile),
	as.character(csvfile),
	as.character(names),
	as.integer(begin),
	as.integer(end),
	ncols, as.integer(verbose), PACKAGE="descr")

    return (invisible(NULL))
}

