# Get data from Opasnet
#
# filename - Name of the file
# wiki - Source Wiki: opasnet_en (default), opasnet_fi, heande (.htaccess protected)
# unzip - File name in package (if compressed)
#
# Returns file contents (loaded using curl)

opasnet.data <- function(filename,wiki='', unzip='') {

	now <- Sys.time()
	
	file <- opasnet.file_url(filename, wiki)
	
	if (unzip != '')
	{
		f <- tempfile(pattern = 'opasnet.data.', fileext = '.zip')
		bin <- getBinaryURL(file)
		con <- file(f, open = "wb")
		writeBin(bin, con)
		close(con)
		con <- unz(f, unzip)
		return(paste(readLines(con),collapse="\n"))
	}
	else
	{	
		return(getURL(file))
	}
}

# Get table data (e.g. csv) from Opasnet
#
# filename - Name of the file
# wiki - Source Wiki: opasnet_en (default), opasnet_fi, heande (.htaccess protected)
# unzip - File name in package (if compressed)
#
# Returns file contents in table (loaded using curl)

opasnet.csv <- function(filename, wiki='', unzip = '', ...) {

	now <- Sys.time()
	
	file <- opasnet.file_url(filename, wiki)

	if (unzip != '')
	{
		f <- tempfile(pattern = 'opasnet.csv.', fileext = '.zip')
		bin <- getBinaryURL(file)
		con <- file(f, open = "wb")
		writeBin(bin, con)
		close(con)
		return(read.table(unz(f, unzip), ...))
	}
	else
	{	
		csv <- getURL(file)
		return(read.table(file = textConnection(csv), ...))
	}
	
	
}

# Get R data from Opasnet
#
# filename - Name of the file
# wiki - Source Wiki: opasnet_en (default), opasnet_fi, heande (.htaccess protected)
# unzip - File name in package (if compressed)
#
# Loads file contents to .GlobalEnv

#opasnet.R <- function(filename,wiki='', unzip='') {
#	
#	now <- Sys.time()
#	
#	file <- opbase.file_url(filename, wiki)
#	
#	if (unzip != '')
#	{
#		f <- tempfile(pattern = 'opasnet.R.', fileext = '.zip')
#		bin <- getBinaryURL(file)
#		con <- file(f, open = "wb")
#		writeBin(bin, con)
#		close(con)
#		con <- unz(f, unzip)
#		load(con, .GlobalEnv)
#		#return(paste(readLines(con),collapse="\n"))
#	}
#	else
#	{
#		load(getURL(file), .GlobalEnv)
#		#return(getURL(file))
#	}
#}

# Private function to get file url for given wiki
opasnet.file_url <- function(filename, wiki)
{
	# Parse arguments
	targs <- strsplit(commandArgs(trailingOnly = TRUE),",")
	args = list()
	if (length(targs) > 0)
		for(i in targs[[1]])
		{
			tmp = strsplit(i,"=")
			key <- tmp[[1]][1]
			value <- tmp[[1]][2]
			args[[key]] <- value
		}
	
	if (wiki == '')
	{
		if (is.null(args$user)) stop('Wiki cannot be resolved!')
		wiki <- args$user
	}
	if (wiki == 'opasnet_en' || wiki == 'op_en')
	{
		file <- paste("http://en.opasnet.org/en-opwiki/images/",filename,sep='')
	}
	if (wiki == 'opasnet_fi' || wiki == 'op_fi')
	{
		file <- paste("http://fi.opasnet.org/fi_wiki/images/",filename,sep='')
	}
	if (wiki == 'heande')
	{
		file <- paste("http://",args$ht_username,":",args$ht_password,"@heande.opasnet.org/heande/images/",filename,sep='')
	}
	return(file)
}


# OPASNET.DATA #####################################
## opasnet.data downloads a file from Finnish Opasnet wiki, English Opasnet wiki, or Opasnet File.
## Parameters: filename is the URL without the first part (see below), wiki is "opasnet_en", "opasnet_fi", or "M-files".
## If table is TRUE then a table file for read.table function is assumed; all other parameters are for this read.table function.
#
#opasnet.data <- function(filename, wiki = "opasnet_en", table = FALSE, ...)
#{
#if (wiki == "opasnet_en") {
#file <- paste("http://en.opasnet.org/en-opwiki/images/", filename, sep = "")
#}
#if (wiki == "opasnet_fi") {
#file <- paste("http://fi.opasnet.org/fi_wiki/images/", filename, sep = "")
#}
#if (wiki == "M-files") {
#file <- paste("http://http://fi.opasnet.org/fi_wiki/extensions/mfiles/", filename, sep = "")
#}
#
#if(table == TRUE) {
#file <- re#ad.table(file, header = FALSE, sep = "", quote = "\"'",
#           dec = ".", row.names, col.names,
#           as.is = !stringsAsFactors,
#           na.strings = "NA", colClasses = NA, nrows = -1,
#           skip = 0, check.names = TRUE, fill = !blank.lines.skip,
#           strip.white = FALSE, blank.lines.skip = TRUE,
#           comment.char = "#",
#           allowEscapes = FALSE, flush = FALSE,
#           stringsAsFactors = default.stringsAsFactors(),
#           fileEncoding = "", encoding = "unknown")
#return(file)
#}
#else {return(ge#tURL(file))}
#}

opasnet.page <- function(pagename, wiki = "") {
	if (wiki == '')
	{
		if (is.null(args$user)) stop('Wiki cannot be resolved!')
		wiki <- args$user
	}
	if (wiki == "opasnet_en" | wiki == "op_en")
	{
		url <- paste("http://en.opasnet.org/en-opwiki/index.php?title=", pagename, sep = "")
	}
	if (wiki == "opasnet_fi" | wiki == "op_fi")
	{
		url <- paste("http://fi.opasnet.org/fi_wiki/index.php?title=", pagename, sep = "")
	}
	if (wiki == 'heande')
	{
		url <- paste("http://",args$ht_username,":",args$ht_password,"@heande.opasnet.org/heande/index.php?title=", pagename, sep = "")
	}
	return(getURL(url))
}