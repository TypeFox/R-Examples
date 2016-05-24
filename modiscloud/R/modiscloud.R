#' @include fermat.R
 
require(sfsmisc)		# for digitsBase
require(raster)			# for raster
#sourcedir = '/Dropbox/_njm/'
#source3 = '_genericR_v1.R'
#source(paste(sourcedir, source3, sep=""))
#roxygenize()


#' Convert to data.frame, without factors
#' 
#' Shortcut for: \code{as.data.frame(x, row.names=NULL, stringsAsFactors=FALSE)}
#' 
#' This function, and \code{\link{adf2}}, are useful for dealing with errors due to 
#' automatic conversion of some columns to factors.  Another solution may be to prepend
#' \code{options(stringsAsFactors = FALSE)} at the start of one's script, to turn off all default stringsAsFactors silliness.
#' 
#' @param x matrix or other object transformable to data.frame
#' @return data.frame
#' @export
#' @seealso \code{\link{adf2}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' adf(x)
adf <- function(x)
	{
	return(as.data.frame(x, row.names=NULL, stringsAsFactors=FALSE))
	}







#' Convert to data.frame, without factors
#' 
#' Shortcut for: \code{tmp_rownames = 1:nrow(x); as.data.frame(x, row.names=tmp_rownames, stringsAsFactors=FALSE)}
#' 
#' This function, and \code{\link{adf2}}, are useful for dealing with errors due to 
#' automatic conversion of some columns to factors.  Another solution may be to prepend
#' \code{options(stringsAsFactors = FALSE)} at the start of one's script, to turn off all default stringsAsFactors silliness.
#'
#' In adf2, rownames are forced to be numbers; this can prevent errors due to e.g. repeated rownames
#' after an \code{rbind} operation.
#'
#' @param x matrix or other object transformable to data.frame
#' @return data.frame
#' @export
#' @seealso \code{\link{adf}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' adf2(x)
adf2 <- function(x)
	{
	# Deals with the problem of repeated row names
	rownames = 1:nrow(x)
	return(as.data.frame(x, row.names=rownames, stringsAsFactors=FALSE))
	}



#######################################################
# unlist_df:
#######################################################
#' Unlist the columns in a data.frame
#' 
#' Sometimes, matrices or data.frames will malfunction due to their having lists as columns
#' and other weirdness.  This is a shortcut for \code{data.frame(lapply(df, function(x) unlist(x)))}.
#' 
#' @param df matrix or other object transformable to data.frame
#' @return data.frame
#' @export
#' @seealso \code{\link{unlist_df2}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' df = adf(matrix(c(1,2,3,4,5,6), nrow=3, ncol=2))
#' df2 = unlist_df(df)
#' df2
#'
unlist_df <- function(df)
	{
	outdf <- data.frame(lapply(df, function(x) unlist(x)))
	}



#######################################################
# unlist_df2:
#######################################################
#' Unlist the columns in a data.frame, with more checks
#' 
#' Sometimes, matrices or data.frames will malfunction due to their having lists as columns
#' and other weirdness. This runs \code{\link{unlist_df}} and additional checks.
#' 
#' @param df matrix or other object transformable to data.frame
#' @return data.frame
#' @export
#' @seealso \code{\link{unlist_df}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' df = adf(matrix(c(1,2,3,4,5,6), nrow=3, ncol=2))
#' df2 = unlist_df2(df)
#' df2
unlist_df2 <- function(df)
	{
	store_colnames = names(df)
	
	outdf = NULL
	
	numrows = dim(df)[1]
	
	for (i in 1:ncol(df))
		{
		#print(names(df)[i])
		tmpcol = unlist(df[, i])
		#print(length(tmpcol))
		
		# Error check; e.g. blank cells might screw it up
		if (length(tmpcol) < numrows)
			{
			tmpcol2 = df[,i]
			tmpcol = as.character(tmpcol2)
			}
		
		outdf = cbind(outdf, tmpcol)
		}
	
	#outdf = adf2(outdf)
	
	names(outdf) = store_colnames
	return(outdf)
	}



#######################################################
# slashslash:
#######################################################
#' Remove double slash (slash a slash)
#' 
#' Shortcut for: \code{gsub(pattern="//", replacement="/", x=tmpstr)}
#' 
#' This function is useful for removing double slashes that can
#' appear in full pathnames due to inconsistencies in trailing
#' slashes in working directories etc.
#'
#' @param tmpstr a path that you want to remove double slashes from
#' @return outstr a string of the fixed path
#' @export
#' @seealso \code{\link{gsub}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' tmpstr = "/Library/Frameworks/R.framework/Versions/2.15/Resources/library/modiscloud/extdata/2002raw//MYD03.A2002185.0645.005.2009192031332.hdf"
#'
#' outstr = slashslash(tmpstr)
#' outstr
#'
slashslash <- function(tmpstr)
	{
	outstr = gsub(pattern="//", replacement="/", x=tmpstr)
	return(outstr)
	}


#######################################################
# byteint2bit:
#######################################################
#' Convert a byte integer (0-255) to a list of 8 bits
#'
#' MODIS HDFs converted to tifs result in each pixel having a 
#' integer value from 0 to 255.  byteint2bit converts this to a
#' string of 8 bits (0/1 values).
#'
#' In the case of MODIS, the bits are read from RIGHT to LEFT, which
#' is unnatural, so by default the function uses rev() to reverse
#' the order of reading.
#'
#' Note: the reversal means that, when interpreting the two 2-bit strings in the
#' MODIS image, the interpretation of 01 and 10 is reversed from in the MODIS 
#' documentation.
#'
#' @param intval, an integer between 0 and 255
#' @param reverse, a logical (TRUE/FALSE) indicating whether or not the string of
#'        bits should be reversed. Default: TRUE
#' @return \code{byte_in_binary}, the 8 bits (0/1 values) in the correct order
#' @export
#' @seealso \code{\link{extract_bit}}
#' @seealso \code{\link{byteint2bit_list}}
#' @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite Ackerman2010
#'   @cite GoldsmithMatzkeDawson2013
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' byteint2bit(intval=0, reverse=FALSE)
#' byteint2bit(intval=0, reverse=TRUE)
#' byteint2bit(intval=1, reverse=FALSE)
#' byteint2bit(intval=1, reverse=TRUE)
#' byteint2bit(intval=254, reverse=FALSE)
#' byteint2bit(intval=254, reverse=TRUE)
#' byteint2bit(intval=255, reverse=FALSE)
#' byteint2bit(intval=255, reverse=TRUE)
#' 
byteint2bit <- function(intval, reverse=TRUE)
	{
	require(sfsmisc)		# for digitsBase

	if (reverse == TRUE)
		{
		byte_in_binary = rev(t(digitsBase(intval, base= 2, 8)))
		} else {
		byte_in_binary = c(t(digitsBase(intval, base= 2, 8)))
		}
	return(byte_in_binary)
	}


#######################################################
# byteint2bit_list: 
#######################################################
#' Convert a list of byte integer (0-255) to a table of 8 bits per row
#'
#' MODIS HDFs converted to tifs result in each pixel having a 
#' integer value from 0 to 255.  byteint2bit converts this to a
#' string of 8 bits (0/1 values).
#'
#' In the case of MODIS, the bits are read from RIGHT to LEFT, which
#' is unnatural, so by default the function uses rev() to reverse
#' the order of reading.
#'
#' This function is just the \code{mapply} (multiply apply) version of \code{\link{byteint2bit}}
#'
#' Note: the reversal means that, when interpreting the two 2-bit strings in the
#' MODIS image, the interpretation of 01 and 10 is reversed from in the MODIS 
#' documentation.
#'
#' @param grdr_vals, a list of integers between 0 and 255
#' @param reverse, a logical (TRUE/FALSE) indicating whether or not the string of
#'        bits should be reversed. Default: TRUE
#' @return \code{byte_in_binary_table}, a table, each cell has the 8 bits (0/1 values) in the correct order
#' @export
#' @seealso \code{\link{extract_bit}}
#' @seealso \code{\link{byteint2bit}}
#' @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite Ackerman2010
#'   @cite GoldsmithMatzkeDawson2013
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' grdr_vals = c(0,1,2,253,254,255)
#' byteint2bit_list(grdr_vals, reverse=TRUE)
#'
byteint2bit_list <- function(grdr_vals, reverse=TRUE)
	{
	if (reverse == TRUE)
		{
		byte_in_binary_table = adf2(t(mapply(FUN=byteint2bit, grdr_vals)))
		} else {
		byte_in_binary_table = adf2(t(mapply(FUN=byteint2bit, grdr_vals, reverse=FALSE)))
		}
		
	# Apply header names for table
	tmpnames = paste("b", seq(0, ncol(byte_in_binary_table)-1, 1), sep="")
	names(byte_in_binary_table) = tmpnames
	
	return(byte_in_binary_table)
	}


#######################################################
# check_for_matching_geolocation_files:
#######################################################
#' Checks that every MODIS cloud project HDF has a matching MOD03 file
#' 
#' Each MOD35_L2 cloud mask product file requires a corresponding
#' MOD03 geolocation file to be successfully processed with the MRTSwath tool.
#'
#' MRTSwath is the MRT (MODIS Reprojection Tool) for the MODIS
#' level 1 and level 2 products (cloud mask is level 2, I think).
#'
#' E.g. this cloud mask file:
#'
#' MOD35_L2.A2012001.0420.005.2012001131638.hdf
#'
#' ...goes with this corresponding geolocation file:
#'
#' MOD03.A2012001.0420.005.2012001104706.hdf
#'
#' ...which is a large file (~30 MB) containing detailed information
#' on the position, tilt, etc. of the MODIS satellite.
#'
#' For whatever reason, even a search done at the same time
#' at http://reverb.echo.nasa.gov/reverb/#utf8=%E2%9C%93&spatial_map=satellite&spatial_type=rectangle
#' will not necessarily return the same number of MOD35_L2 and MOD03 files. 
#' MRTSwath tool needs one of each, however.
#'
#' @param moddir the string describing the directory containing the MOD35_L2 and MOD03 files; both must be in the same directory.  Default: getwd(), which gives the present working directory.
#' @param modtxt the text string indicating which HDF files are the MODIS cloud product (or hypothetically, other product). Default: MOD35_L2 (MODIS cloud mask product)
#' @param geoloctxt the text string indicating which HDF files are the MODIS geolocation files (or hypothetically, another set of files). Default: MOD03
#' @param return_geoloc if TRUE, return the list of unmatched geolocation files (e.g. MOD03 or MYD03)
#' @param return_product if TRUE, return the list of unmatched product files (e.g. MOD35_L2 or MYD35_L2 or MOD06_L2 or MYD06_L2)
#' @return data.frame of matching files; or a list of non-matching files, if \code{return_geoloc} or \code{return_product} are TRUE.
#' @export
#' @seealso \code{\link{extract_time_from_MODISfn}}
#' @seealso \code{\link{modfns_to_ftp_download_cmd}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # Check your working directory
#' moddir = getwd()
#' 
#' # Here are some example MODIS files in modiscloud/extdata/
#' # Code excluded from CRAN check because it depends on modiscdata
#' \dontrun{
#' library(devtools)
#' # The modiscdata (MODIS c=cloud data=data) package is too big for CRAN (60 MB); so it is available on github:
#' # https://github.com/nmatzke/modiscdata
#' # If we can't get install_github() to work, try install_url():
#' # install_github(repo="modiscdata", username="nnmatzke")
#' install_url(url="https://github.com/nmatzke/modiscdata/archive/master.zip")
#' library(modiscdata)
#' moddir = system.file("extdata/2002raw/", package="modiscdata")
#'
#' # You need to have some e.g. MOD files in it (from the MODIS-TERRA platform)
#' # (*won't* work with the default files stored in modiscloud/extdata/2002raw/)
#' list.files(path=moddir, pattern="MOD")
#'
#' # This directory actually has MYD files (from the MODIS-AQUA platform)
#' # (*will* work with the default files stored in modiscloud/extdata/2002raw/)
#' list.files(path=moddir, pattern="MYD")
#'
#' # Check for matches (for MODIS-TERRA platform)
#' # (*won't* work with the default files stored in modiscloud/extdata/2002raw/)
#' check_for_matching_geolocation_files(moddir=moddir, modtxt="MOD35_L2", geoloctxt="MOD03", return_geoloc=FALSE, return_product=FALSE)
#'
#' # Check for matches (for MODIS-AQUA platform)
#' # (*will* work with the default files stored in modiscloud/extdata/2002raw/)
#' check_for_matching_geolocation_files(moddir=moddir, modtxt="MYD35_L2", geoloctxt="MYD03", return_geoloc=FALSE, return_product=FALSE)
#' }
#'
check_for_matching_geolocation_files <- function(moddir=getwd(), modtxt="MOD35_L2", geoloctxt="MOD03", return_geoloc=FALSE, return_product=FALSE)
	{

	# Get fns with suffix, from either directory or fns list
	# Do NOT use dot in the suffix
	get_fns_matching_txt <- function(tmpdir=NULL, fns=NULL, text_to_match = NULL, returnfullnames=TRUE)
		{
		# If tmpdir is NOT null, get those files from list.files.
		if (is.null(tmpdir) == FALSE)
			{
			fns = list.files(tmpdir, full.names=returnfullnames)
			fns_without_paths = list.files(tmpdir, full.names=FALSE)
			}
		
		# Return true/false for matched text
		TF = grepl(pattern=text_to_match, x=fns_without_paths)
		
		matching_fns = fns[TF]
		return(matching_fns)
		}



	mod03_fns = sort(slashslash(get_fns_matching_txt(moddir, text_to_match=geoloctxt)))
	mod35_L2_fns = sort(slashslash(get_fns_matching_txt(moddir, text_to_match=modtxt)))
	#head(mod03_fns)
	#head(mod35_L2_fns)

	# Check if you have the right # of files
	if (length(mod03_fns) == length(mod35_L2_fns))
		{
		cat("\nProcessing ", length(mod03_fns), " files.\n", sep="")
		
		# Return the list, sorted
		fns_df = cbind(mod35_L2_fns, mod03_fns)
		fns_df = adf(fns_df)
		
		return(fns_df)
		
		} else {
		cat("\nWARNING: You have ", length(mod35_L2_fns), " ", modtxt, " files & ", length(mod03_fns), " ", geoloctxt, " files.\nWill attempt to find just the matching files", sep="")
		}
	
	
	# Get the datestring for each MOD03 file
	mod03_idstrings = rep(NA, times = length(mod03_fns))
	for (i in 1:length(mod03_fns))
		{
		words = strsplit(x=mod03_fns[i], split="\\.")[[1]]
		idnums = words[2:4]
		mod03_idstrings[i] = paste(idnums, sep="", collapse=".")
		
		}
	#mod03_idstrings

	# Get the datestring for each mod35_L2 file
	mod35_L2_idstrings = rep(NA, times = length(mod35_L2_fns))
	for (i in 1:length(mod35_L2_fns))
		{
		words = strsplit(x=mod35_L2_fns[i], split="\\.")[[1]]
		idnums = words[2:4]
		mod35_L2_idstrings[i] = paste(idnums, sep="", collapse=".")
		}
	#mod35_L2_idstrings
	
	
	# Find which match
	mod35_L2_in_mod03_TF = mod35_L2_idstrings %in% mod03_idstrings
	mod35_L2_fns_dropped = mod35_L2_fns[mod35_L2_in_mod03_TF == FALSE]
	mod35_L2_fns = mod35_L2_fns[mod35_L2_in_mod03_TF == TRUE]
	mod35_L2_idstrings = mod35_L2_idstrings[mod35_L2_in_mod03_TF == TRUE]
	
	mod03_in_mod35_L2_TF = mod03_idstrings %in% mod35_L2_idstrings
	mod03_fns_dropped = mod03_fns[mod03_in_mod35_L2_TF == FALSE]
	mod03_fns = mod03_fns[mod03_in_mod35_L2_TF == TRUE]
	mod03_idstrings = mod03_idstrings[mod03_in_mod35_L2_TF == TRUE]
	
	
	# Correct // to / (if any)
	# Could also use slashslash
	mod03_fns = gsub(pattern="//", replacement="/", x=mod03_fns)
	mod35_L2_fns = gsub(pattern="//", replacement="/", x=mod35_L2_fns)
	mod35_L2_fns_dropped = gsub(pattern="//", replacement="/", x=mod35_L2_fns_dropped)
	mod03_fns_dropped = gsub(pattern="//", replacement="/", x=mod03_fns_dropped)
	
	# Check lengths (manual)
	length(mod35_L2_fns)
	length(mod03_fns)	
	length(mod35_L2_idstrings)
	length(mod03_idstrings)
	sum(mod35_L2_idstrings == mod03_idstrings)
	
	# Return the list or matching files, sorted
	fns_df = cbind(mod35_L2_fns, mod03_fns)
	fns_df = adf(fns_df)
	
	# Print the dropped files
	cat("\n", sep="")
	cat("\nWarning: ", length(mod35_L2_fns_dropped), " ", modtxt, " files dropped with no matches:\n", sep="")
	cat(head(mod35_L2_fns_dropped), sep="\n")
	cat("...", sep="")
	cat("\n", sep="")

	cat("\nWarning: ", length(mod03_fns_dropped), " ", geoloctxt, " files dropped with no matches:\n", sep="")
	cat(head(mod03_fns_dropped), sep="\n")
	cat("...", sep="")
	cat("\n", sep="")
	
	# Return unmatched geolocation files, if desired
	if (return_geoloc == TRUE)
		{
		return(mod03_fns_dropped)
		}

	# Return unmatched product files, if desired
	if (return_product == TRUE)
		{
		return(mod35_L2_fns_dropped)
		}
	
	# Otherwise, return the matched files...
	return(fns_df)	
	}






#######################################################
# dates_from_fileslist:
#######################################################
#' Convert each MODIS filename to year, month, day, hourmin
#'
#' For each MODIS filename in a list, return a table of date/time information.
#'
#' @param fns filenames
#' @return \code{dates_table} a table of dates
#' @export
#' @seealso \code{\link{yearday_to_date}}
#' @seealso \code{\link{extract_time_from_MODISfn}}
#' @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite Ackerman2010
#'   @cite GoldsmithMatzkeDawson2013
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' fns = c("MOD03.A2008001.0400.005.2010216170200.hdf", "MOD03.A2007002.0300.005.2010216170200.hdf")
#' dates_from_fileslist(fns)
#'
dates_from_fileslist <- function(fns)
	{
	require(date)
	
	# Get the year, day, hourmin
	yd_table = adf2(t(mapply(FUN=extract_time_from_MODISfn, fns)))
	
	# Get the month
	ym_table = adf2(t(mapply(FUN=yearday_to_date, yd_table$year, yd_table$day)))
	
	dates_table = cbind(ym_table$year, ym_table$month, ym_table$day, yd_table$day, yd_table$hourmin)
	
	dates_table = unlist_df2(dates_table)
	
	dates_table = adf2(dates_table)
	names(dates_table) = c("year", "month", "day", "julian", "hourmin")
	
	# Convert to numeric
	dates_table$year = as.numeric(dates_table$year)
	dates_table$month = as.numeric(dates_table$month)
	dates_table$day = as.numeric(dates_table$day)
	dates_table$julian = as.numeric(dates_table$julian)

	POSIX_ct_date = mapply(make_POSIXct_date, year=dates_table$year, month=dates_table$month, day=dates_table$day, hourmin=dates_table$hourmin)
	
	dates_table = cbind(dates_table, POSIX_ct_date)
	
	#head(dates_table)
	#tail(dates_table)
	
	return(dates_table)
	}





#######################################################
# extract_bit:
#######################################################
#' Get the value of a particular bit in a byte
#'
#' In many MODIS products, the information for each pixel is stored in a byte, which 
#' represents numbers between 0 and 255.  A byte is made up of 8 bits, and is a 
#' base 2 representation of the number.  In the MODIS cloud product, each bit encodes
#' information; see documentation of the MODIS cloud product in question.
#' 
#' This function takes a byte and extracts the bits.  The bits are in whatever order they
#' are in the byte.  If reading in the reverse direction is needed, the user should use \code{\link{rev}}.
#' 
#' @param intval An integer between 0 and 255
#' @param bitnum The bit number to return, between 1 and 8
#' @return \code{bitval} The value of the bit (0 or 1). 
#' @export
#' @seealso \code{\link{byteint2bit}}
#' @seealso \code{\link{get_bitgrid_2bits}}
#' @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite Ackerman2010
#'   @cite GoldsmithMatzkeDawson2013
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' intval = 123
#' extract_bit(intval, bitnum=1)
#' extract_bit(intval, bitnum=2)
#' extract_bit(intval, bitnum=8)
#'
extract_bit <- function(intval, bitnum)
	{
	require(sfsmisc)		# for digitsBase
	
	# Convert to binary, read correctly (MODIS reads right-to-left)
	binval = rev(t(digitsBase(intval, base= 2, 8)))
	
	bitval = binval[bitnum]
	return(bitval)
	}



#######################################################
# extract_fn_from_path:
#######################################################
#' Get the filename from a path
#'
#' The filename is split on slashes, and the last item is taken; this should be just
#' the filename.
#' 
#' @param fn_with_path The filename, with partial or full path
#' @return \code{fn} The extracted filename
#' @export
#' @seealso \code{\link[base]{strsplit}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' fn_with_path = "/Library/Frameworks/R.framework/Versions/2.15/Resources/library/modiscloud/extdata/2002raw/MYD35_L2.A2002185.1910.005.2007206043609.hdf"
#' extract_fn_from_path(fn_with_path)
#'
extract_fn_from_path <- function(fn_with_path)
	{
	words = strsplit(fn_with_path, split="/")[[1]]
	fn = words[length(words)]
	return(fn)
	}




#######################################################
# extract_time_from_MODISfn:
#######################################################
#'  Extract the year, day, and hour from a MODIS filename
#'
#' The filenames of MODIS images contain the following date
#' information: MOD35_L2.Ayyyyddd.hhhh.etc.
#'
#'     MODLAND Level 2 products\cr
#'     ESDT.AYYYYDDD.HHMM.CCC.YYYYDDDHHMMSS.hdf\cr
#'     
#'     ESDT = Earth Science Data Type name (e.g., MOD14)\cr
#'     YYYYDDD = MODIS acquisition year and Julian day\cr
#'     HHMM = MODIS acquisition UTC time\cr
#'     CCC = Collection number\cr
#'     YYYYDDDHHMMSS = Processing Year, Julian day and UTC Time\cr
#'     hdf = Suffix denoting HDF file\cr
#'
#' DDD is the day of the year, from 001 to 365 (or 366 for leap years)
#' 
#' 
#' @param fn The MODIS filename of interest
#' @return \code{newdate} A list with three items: $month, $day, $year
#' @export
#' @seealso \code{\link{dates_from_fileslist}}
#' @seealso \url{http://landweb.nascom.nasa.gov/cgi-bin/QA_WWW/newPage.cgi?fileName=hdf_filename}
#'   @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite NASA2001
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' yearday_to_date(year=2012, day=364)
#' # $month
#' # [1] 12
#' # 
#' # $day
#' # [1] 29
#' # 
#' # $year
#' # [1] 2008
#'
#' fn = "MOD03.A2008001.0400.005.2010216170200.hdf"
#' extract_time_from_MODISfn(fn)
#'
extract_time_from_MODISfn <- function(fn)
	{
	#tmpstr = "MOD03.A2008001.0400.005.2010216170200.hdf"
	tmpstr = fn
	
	# Split filename on "."
	words = strsplit(tmpstr, split="\\.")
	
	yearday = words[[1]][2]
	year = as.numeric(substr(yearday, start=2, stop=5))
	day = as.numeric(substr(yearday, start=6, stop=8))
	
	# Keep as txt for now
	hourmin = words[[1]][3]
	
	datelist = NULL
	datelist$year = year
	datelist$day = day
	datelist$hourmin = hourmin
	
	return(datelist)
	}



#' Test function
#'
#' This is an example function for making R packages and R documentation
#' with \code{\link[roxygen2:roxygenize]{roxygen2}} and \code{\link[roxygen2:roxygenize]{roxygenize}}.
#' 
#' This function just calculates the mean.
#' 
#' @param d hi
#' @return meand
#' @export
#' @seealso \code{\link{fermat.test}}
#' @references refline1
#' refline2
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' foo(c(1,2,3,4)) # 2.5
#'
foo <- function(d)
	{
	meand = mean(d)
	return(meand)
	}



#' Test function
#'
#' This is an example function for making R packages and R documentation
#' with \code{\link[roxygen2:roxygenize]{roxygen2}} and \code{\link[roxygen2:roxygenize]{roxygenize}}.
#' 
#' This function just calculates the mean.
#' 
#' Summary_paragraph
#' 
#' Details_paragraph
#'
#' @param d description_of_input_param
#' @return description_of_what_is_returned
#' @export
#' @seealso \code{\link{fermat.test}}
#' @references refline1
#' refline2
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' foo(c(1,2,3,4)) # 2.5
#'
foo2 <- function(d)
	{
	meand = mean(d)
	return(meand)
	}








#######################################################
# get_bitgrid:
#######################################################
#' Extract the value of a particular bit at each pixel
#'
#' Given an input grid (as a SpatialGridDataFrame object), extract a single desired bit value from
#' all of the pixels
#'
#' Note: this means that, when interpreting the two 2-bit strings in the
#' MODIS image, you would have to extract each bit separately.
#'
#' @param grd A SpatialGridDataFrame, derived from e.g. readGDAL
#' @param bitnum The bit you would like to extract
#' @return \code{grdr_vals_bits}, a list of 0/1 values for the bit in question (numeric)
#' @export
#' @seealso \code{\link{extract_bit}}
#' @seealso \code{\link{get_bitgrid_2bits}}
#' @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite Ackerman2010
#'   @cite GoldsmithMatzkeDawson2013
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' #######################################################
#' # Load a TIF
#' #######################################################
#' # Code excluded from CRAN check because it depends on modiscdata
#' \dontrun{
#' library(devtools)
#' # The modiscdata (MODIS c=cloud data=data) package is too big for CRAN (60 MB); so it is available on github:
#' # https://github.com/nmatzke/modiscdata
#' # If we can't get install_github() to work, try install_url():
#' # install_github(repo="modiscdata", username="nnmatzke")
#' install_url(url="https://github.com/nmatzke/modiscdata/archive/master.zip")
#' library(modiscdata)
#' tifsdir = system.file("extdata/2002tif/", package="modiscdata")
#' tiffns = list.files(tifsdir, pattern=".tif", full.names=TRUE)
#' tiffns
#' 
#' library(rgdal)	# for readGDAL
#' 
#' # numpixels in subset
#' xdim = 538
#' ydim = 538
#' 
#' 
#' # Read the grid and the grid metadata		
#' coarsen_amount = 1
#' xdim_new = xdim / floor(coarsen_amount)
#' ydim_new = ydim / floor(coarsen_amount)
#' 		
#' fn = tiffns[1]
#' grd = readGDAL(fn, output.dim=c(ydim_new, xdim_new))
#' 
#' grdproj = CRS(proj4string(grd))
#' grdproj
#' grdbbox = attr(grd, "bbox")
#' grdbbox
#' 
#' #######################################################
#' # Extract a particular bit for all the pixels in the grid
#' #######################################################
#' bitnum = 2
#' grdr_vals_bits = get_bitgrid(grd, bitnum)
#' length(grdr_vals_bits)
#' grdr_vals_bits[1:50]
#' }
#' 
get_bitgrid <- function(grd, bitnum)
	{
	# Convert to raster, and extract grid values	
	grdr = raster(grd)
	pointscount_on_SGDF_points = coordinates(grd)
	grdr_vals = extract(grdr, pointscount_on_SGDF_points)
	
	# Get unique values on the grid, and convert each to bits
	uniq_vals = sort(unique(grdr_vals))
	uniq_vals_bits = mapply(FUN=extract_bit, intval=uniq_vals, bitnum=1)

	# Replace the values
	grdr_vals_bits = grdr_vals
	for (j in 1:length(uniq_vals))
		{
		tmpval = uniq_vals[j]
		newval = uniq_vals_bits[j]
		
		grdr_vals_bits[grdr_vals == tmpval] = newval
		}

	return(grdr_vals_bits)
	}

#######################################################
# get_bitgrid_2bits: 
#######################################################
#' extract the value of 2 particular bits at each pixel
#'
#' Some MODIS cloud products have 2-bit codes, e.g. the MOD35_L2 product has a 2-bit
#' code specifying the 4 possible outputs of the cloud-classification algorithm:
#' confident cloudy, probable cloudy, probable clear, confident clear.
#'
#' Note: this means that, when interpreting the two 2-bit strings in the
#' MODIS image, you have to reverse the bit order. This is done by default here.
#'
#' @param grd A SpatialGridDataFrame, derived from e.g. readGDAL
#' @param bitnums The bit you would like to extract
#' @return \code{grdr_vals_bitstrings}, a list of strings for the bits in question (string, e.g. "11" or "01" or "00")
#' @export
#' @seealso \code{\link{extract_bit}}
#' @seealso \code{\link{get_bitgrid}}
#' @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite Ackerman2010
#'   @cite GoldsmithMatzkeDawson2013
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' #######################################################
#' # Load a TIF
#' #######################################################
#' # Code excluded from CRAN check because it depends on modiscdata
#' \dontrun{
#' library(devtools)
#' # The modiscdata (MODIS c=cloud data=data) package is too big for CRAN (60 MB); so it is available on github:
#' # https://github.com/nmatzke/modiscdata
#' # If we can't get install_github() to work, try install_url():
#' # install_github(repo="modiscdata", username="nnmatzke")
#' install_url(url="https://github.com/nmatzke/modiscdata/archive/master.zip")
#' library(modiscdata)
#' tifsdir = system.file("extdata/2002tif/", package="modiscdata")
#' tiffns = list.files(tifsdir, pattern=".tif", full.names=TRUE)
#' tiffns
#' 
#' library(rgdal)	# for readGDAL
#' 
#' # numpixels in subset
#' xdim = 538
#' ydim = 538
#' 
#' 
#' # Read the grid and the grid metadata		
#' coarsen_amount = 1
#' xdim_new = xdim / floor(coarsen_amount)
#' ydim_new = ydim / floor(coarsen_amount)
#' 		
#' fn = tiffns[1]
#' grd = readGDAL(fn, output.dim=c(ydim_new, xdim_new))
#' 
#' grdproj = CRS(proj4string(grd))
#' grdproj
#' grdbbox = attr(grd, "bbox")
#' grdbbox
#' 
#' #######################################################
#' # Extract a particular pair of bits for all the pixels in the grid
#' #######################################################
#' bitnum = 2
#' grdr_vals_bitstrings = get_bitgrid_2bits(grd, bitnum)
#' length(grdr_vals_bitstrings)
#' grdr_vals_bitstrings[1:50]
#' }
#' 
get_bitgrid_2bits <- function(grd, bitnums)
	{
	# Convert to raster, and extract grid values	
	grdr = raster(grd)
	pointscount_on_SGDF_points = coordinates(grd)
	grdr_vals = extract(grdr, pointscount_on_SGDF_points)

	# Get unique byte / integer values
	uniq_vals = sort(unique(grdr_vals))
	

	# Get the first bit
	# Replace the values
	grdr_vals_bitA = grdr_vals

	# Get unique values on the grid, and convert each to bits
	uniq_vals_bits = mapply(FUN=extract_bit, intval=uniq_vals, bitnum=bitnums[1])

	for (j in 1:length(uniq_vals))
		{
		tmpval = uniq_vals[j]
		newval = uniq_vals_bits[j]
		
		grdr_vals_bitA[grdr_vals == tmpval] = newval
		}

	# Get the second bit
	# Replace the values
	grdr_vals_bitB = grdr_vals

	# Get unique values on the grid, and convert each to bits
	uniq_vals_bits = mapply(FUN=extract_bit, intval=uniq_vals, bitnum=bitnums[2])

	for (j in 1:length(uniq_vals))
		{
		tmpval = uniq_vals[j]
		newval = uniq_vals_bits[j]
		
		grdr_vals_bitB[grdr_vals == tmpval] = newval
		}

	grdr_vals_bitstrings = paste(grdr_vals_bitA, grdr_vals_bitB, sep="")

	return(grdr_vals_bitstrings)
	}




#######################################################
# get_date_from_POSIXct: 
#######################################################
#' Get the time information from a POSIXct time
#'
#' This function gets the year, month, day, and hourmin, from a POSIXct date (which is the # of seconds from "1970-01-01 00:00.00 UTC").  Returns a formatted date in data.frame format
#'
#' The function contains a check for multiple dates.
#' 
#' @param POSIX_ct_date (# of seconds from "1970-01-01 00:00.00 UTC")
#' @return tdf a time-data.frame containing year, month, day, julian day, and hourmin / sec
#' @export
#' @seealso \code{\link{get_dates_from_POSIXct}}
#' @seealso \code{\link{make_POSIXct_date}}
#' @seealso \code{\link[base]{as.POSIXct}}
#' @seealso \code{\link[base]{as.POSIXlt}}
#' @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite Ackerman2010
#'   @cite GoldsmithMatzkeDawson2013
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # 10 million seconds from 1/1/1970
#' POSIX_ct_date = 10000000
#' get_date_from_POSIXct(POSIX_ct_date)
#' 
get_date_from_POSIXct <- function(POSIX_ct_date)
	{
	# If more than 1 date
	if (length(POSIX_ct_date) > 1)
		{
		POSIX_ct_dates = POSIX_ct_date
		return(get_dates_from_POSIXct(POSIX_ct_dates))
		}	
	
	POSIX_lt_date = as.POSIXlt(POSIX_ct_date, tz="UTC", origin="1970-01-01 00:00.00 UTC")
	
	# If just 1 date
	words = strsplit(x=as.character(POSIX_lt_date), split=" ")[[1]]
	date_char = words[1]
	hourminsec_char = words[2]
	#tz = words[3]
	
	# Extract date
	datenums = as.numeric(strsplit(date_char, split="-")[[1]])
	year = datenums[1]
	month = datenums[2]
	day = datenums[3]
	
	# Extract time
	hourminsec_words = strsplit(hourminsec_char, split=":")[[1]]
	hourchar = hourminsec_words[1]
	minchar = hourminsec_words[2]
	secchar = hourminsec_words[3]
	
	hournum = as.numeric(hourchar)
	minnum = as.numeric(minchar)
	secnum = as.numeric(secchar)
	
	# Get hourmin string
	hourmin = paste(hourchar, minchar, sep="")
	
	# Get Julian day (day count from Jan. 1)
	julian_day = yearmonthday_to_julianday(year, month, day)

	# Make output table
	dtf = cbind(year, month, day, julian_day, hourmin, hournum, minnum, secnum)
	dtf = adf2(dtf)
	names(dtf) = c("year", "month", "day", "julian_day", "hourmin", "hournum", "minnum", "secnum")
	
	# Force all columns to numeric, except for hourmin
	#cls.df(dtf)
	dtf[,-5] = as.numeric(dtf[,-5])
	#cls.df(dtf)
	
	return(dtf)
	}



#######################################################
# get_dates_from_POSIXct:
#######################################################
#' Get the time information from POSIXct times
#'
#' This function, for each input date, gets the year, month, day, and hourmin from a POSIXct date
#' (which is the # of seconds from "1970-01-01 00:00.00 UTC").  Returns a formatted date in data.frame format
#'
#' The function contains a check for singe dates. This is a version of get_dates_from_POSIXct for
#' use with a list of input POSIXct dates,
#' rather than a single date.
#'
#' @param POSIX_ct_dates list of POSIXct dates
#' @return \code{dtf} data.frame containing year, month, day, julian day, and hourmin / sec
#' @export
#' @seealso \code{\link{get_date_from_POSIXct}}
#' @seealso \code{\link{make_POSIXct_date}}
#' @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite Ackerman2010
#'   @cite GoldsmithMatzkeDawson2013
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' 
#' POSIX_ct_dates = c(10000000, 20000000)
#' get_dates_from_POSIXct(POSIX_ct_dates)
#' 
get_dates_from_POSIXct <- function(POSIX_ct_dates)
	{
	# If just 1 date
	if (length(POSIX_ct_dates) == 1)
		{
		POSIX_ct_date = POSIX_ct_dates
		return(get_date_from_POSIXct(POSIX_ct_date))
		}	
	
	POSIX_lt_dates = as.POSIXlt(POSIX_ct_dates, tz="UTC", origin="1970-01-01 00:00.00 UTC")
	
	# If more than one date
	words = unlist(strsplit(x=as.character(POSIX_lt_dates), split=" "))
	words_matrix = matrix(words, ncol=2, byrow=TRUE)
	
	date_char = words_matrix[,1]
	hourminsec_char = words_matrix[,2]
	#tz = words[3]
	
	# Extract date
	datenums = matrix(unlist(strsplit(date_char, split="-")), ncol=3, byrow=TRUE)
	year = as.numeric(datenums[,1])
	month = as.numeric(datenums[,2])
	day = as.numeric(datenums[,3])
	
	# Extract time
	hourminsec_words = matrix(unlist(strsplit(hourminsec_char, split=":")), ncol=3, byrow=TRUE)
	hourchar = hourminsec_words[,1]
	minchar = hourminsec_words[,2]
	secchar = hourminsec_words[,3]
	
	hournum = as.numeric(hourchar)
	minnum = as.numeric(minchar)
	secnum = as.numeric(secchar)
	
	# Get hourmin string
	hourmin = paste(hourchar, minchar, sep="")
	
	# Get Julian day (day count from Jan. 1)
	julian_day = mapply(FUN=yearmonthday_to_julianday, year, month, day)
	
	
	# THIS ORDER OF OPERATIONS TO GET A DATAFRAME NOT A MATRIX!!
	
	# Make output table
	dtf = cbind(year, month, day, julian_day, hourmin, hournum, minnum, secnum)
	
	# unlist the dtf?
	dtf = unlist_df2(dtf)
	dtf = adf2(dtf)

	names(dtf) = c("year", "month", "day", "julian_day", "hourmin", "hournum", "minnum", "secnum")
	
	# Force all columns to numeric, except for hourmin
	#cls.df(dtf)
	dtf$year = as.numeric(dtf$year)
	dtf$month = as.numeric(dtf$month)
	dtf$day = as.numeric(dtf$day)
	dtf$julian_day = as.numeric(dtf$julian_day)
	dtf$hournum = as.numeric(dtf$hournum)
	dtf$minnum = as.numeric(dtf$minnum)
	dtf$secnum = as.numeric(dtf$secnum)
	#cls.df(dtf)
	
	
	return(dtf)
	}



#' Check an integer for pseudo-primality to an arbitrary precision
#' 
#' A number is pseudo-prime if it is probably prime, the basis
#' of which is the probabilistic Fermat test; if it passes two
#' such tests, the chances are better than 3 out of 4 that
#' \eqn{n} is prime.
#'
#' This is an example function for making R packages and R documentation
#' with \code{\link[roxygen2:roxygenize]{roxygen2}} and \code{\link[roxygen2:roxygenize]{roxygenize}}.
#' 
#' @param n the integer to test for pseudoprimality.
#' @param times the number of Fermat tests to perform
#' @return Whether the number is pseudoprime
#' @export
#' @seealso \code{\link{fermat.test}}
#' @references Abelson, Hal; Jerry Sussman, and Julie Sussman.
#' Structure and Interpretation of Computer Programs.
#' Cambridge: MIT Press, 1984.
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' is.pseudoprime(13, 4)# TRUE most of the time
#' 
is.pseudoprime <- function(n, times)
	{
	if (times == 0) TRUE
	else if (fermat.test(n)) is.pseudoprime(n, times - 1)
	else FALSE
	}


#' Check an integer for pseudo-primality to an arbitrary
#' precision, second version
#'
#' A number is pseudo-prime if it is probably prime, the basis
#' of which is the probabilistic Fermat test; if it passes two
#' such tests, the chances are better than 3 out of 4 that
#' \eqn{n} is prime.
#'
#' This is an example function for making R packages and R documentation
#' with \code{\link[roxygen2:roxygenize]{roxygen2}} and \code{\link[roxygen2:roxygenize]{roxygenize}}.
#'
#' @param n the integer to test for pseudoprimality.
#' @param times the number of Fermat tests to perform
#' @return Whether the number is pseudoprime
#' @export
#' @seealso \code{\link{fermat.test}}
#' @references Abelson, Hal; Jerry Sussman, and Julie Sussman.
#' Structure and Interpretation of Computer Programs.
#' Cambridge: MIT Press, 1984.
#' @author Peter Danenberg \email{pcd@@roxygen.org}
#' @examples
#' is.pseudoprime2(13, 4)# TRUE most of the time
#'
is.pseudoprime2 <- function(n, times)
	{
	if (times == 0) TRUE
	else if (fermat.test(n)) is.pseudoprime2(n, times - 1)
	else FALSE
	}




#######################################################
# make_POSIXct_date:
#######################################################
#' Take year, month, day, hourmin, convert to POSIXct
#'
#' Make a standard-formatted date.
#'
#' @param year calendar year
#' @param month number of the month (1-12)
#' @param day number of the day (1-31)
#' @param hourmin a string (derived from MODIS filename) with HHSS (hours, seconds)
#' @param timezone_txt the timezone; default 'UTC', your results may vary for the other timezones R knows about
#' @return \code{POSIX_ct_date}, a POSIXct-formatted date
#' @export
#' @seealso \code{\link{get_date_from_POSIXct}}
#' @seealso \code{\link{get_dates_from_POSIXct}}
#' @seealso \code{\link[base]{as.POSIXct}}
#' @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite Ackerman2010
#'   @cite GoldsmithMatzkeDawson2013
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # Example use of make_POSIXct_date
#' year=2008
#' month=3
#' day=3
#' hourmin="0404"
#' make_POSIXct_date(year, month, day, hourmin, timezone_txt="UTC")
#'
make_POSIXct_date <- function(year, month, day, hourmin, timezone_txt="UTC")
	{
	test = '
	year=2008
	month=3
	day=3
	hourmin="0404"
	'

	month = sprintf("%02.0f", month)
	day = sprintf("%02.0f", day)
	
	hour = substr(hourmin, start=1, stop=2)
	minute = substr(hourmin, start=3, stop=4)
	sec = "00"

	datestr = paste(year, month, day, sep="-")
	timestr = paste(hour, minute, sec, sep=":")
	
	fullstr = paste(datestr, timestr, "UTC", sep=" ")
	
	# Require that the date be interpreted as a UTC date!
	POSIX_ct_date = as.POSIXct(strptime(fullstr, '%Y-%m-%d %H:%M:%S'), tz=timezone_txt)
	
	return(POSIX_ct_date)
	}






#######################################################
# make_cloudcount_product: 
#######################################################
#' Take a cloudcount raster and a number-of-observations raster and make a fraction cloud product
#'
#' Takes a cloudcount raster, and a number-of-observations raster, and makes a fraction cloud product.
#' 
#' @param master_bitgridvalsCC Cloud count raster
#' @param master_bitgridvals0 Raster indicating the count of valid observations
#' @param num_its The total number of images or files, i.e. the number of attempts to see clouds
#' @param grd An input grid with the appropriate projection and extent to match the rasters
#' @param xdim_new The x-dimension (in # of pixels; i.e. # of columns) of the input/output grids
#' @param ydim_new The y-dimension (in # of pixels; i.e. # of rows) of the input/output grids
#' @param countzeros Should zeros be counted? TRUE or FALSE
#' @return \code{grd_final}, a grid with fraction cloudy
#' @export
#' @seealso \code{\link{sum_bitgrid}}
#' @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite Ackerman2010
#'   @cite GoldsmithMatzkeDawson2013
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # For examples, see the function numslist_to_grd()
#' 
make_cloudcount_product <- function(master_bitgridvalsCC, master_bitgridvals0, num_its, grd, xdim_new, ydim_new, countzeros=FALSE)
	{
	require(raster)
	
	# Get the grd projection
	grdproj = CRS(proj4string(grd))
	
	# Subtract the zeros that are due to bit0=0
	numdays_w_obs = master_bitgridvals0
	numdays_wo_obs = num_its-master_bitgridvals0
	
	if (countzeros == TRUE)
		{
		cloudcount = num_its-master_bitgridvalsCC
		} else {
		cloudcount = master_bitgridvalsCC
		}	
	
	cloudcount = cloudcount - numdays_wo_obs
	cloud_percent = cloudcount / numdays_w_obs
	
	unique(numdays_wo_obs)
	unique(cloudcount)
	unique(cloud_percent)
	#hist(cloud_percent)
	
	
	tmpext2 = extent(grd)
	numcells_long = xdim_new
	numcells_lat = ydim_new
	
	xmin = attr(tmpext2, "xmin")
	xmax = attr(tmpext2, "xmax")
	ymin = attr(tmpext2, "ymin")
	ymax = attr(tmpext2, "ymax")
	xcoords = seq(from=xmin, to=xmax, by=((xmax-xmin)/numcells_long))
	ycoords = seq(from=ymax, to=ymin, by=(-1*(ymax-ymin)/numcells_lat))
	xcoords = xcoords[1:(length(xcoords)-1)]
	ycoords = ycoords[1:(length(ycoords)-1)]



	
	# Directions from here:
	# http://r-sig-geo.2731867.n2.nabble.com/Need-fast-method-to-find-indices-of-cells-in-a-3-D-grid-that-contain-x-y-z-coordinates-td2763738.html
	
	# Number the cells in the 3-D array.
	numcols = length(xcoords)
	numrows = length(ycoords)
	#cell_ids <- array(1:(numcols*numrows), dim = c(numrows,numcols))
	cell_ids <- array(cloud_percent, dim = c(numrows,numcols))
	dim(cell_ids)
	#cell_xs <- seq(from=xmin, to=xmax, by=((xmax-xmin)/numcells_long))
	#cell_ys <- seq(from=ymax, to=ymin, by=((ymax-ymin)/numcells_long))
	
	# Get the numbers of the grid cells containing the coordinates x, y, z.
	# x, y, z are the coordinates I want to look up.  
	# xnb, ynb, and znb are the boundaries of the cells in the x, y, and z
	# directions, respectively.  
	#col_ind <- findInterval(xy_points$x, xcoords, all.inside = T)
	#row_ind <- findInterval(xy_points$y, ycoords, all.inside = T)
	#index.array <- array(c(row_ind, col_ind), dim=c(length(xy_points$x), 2))
	#indices_you_want = cell_ids[index.array]  # the numbers of the cells I want 
	
	
	xcoords1 = matrix(data=xcoords, nrow=numrows, ncol=numcols, byrow=FALSE)
	ycoords1 = matrix(data=ycoords, nrow=numrows, ncol=numcols, byrow=TRUE)
	#xcoords1
	#ycoords1
	
	xy_raster = NULL
	xy_raster$x = c(xcoords1)
	xy_raster$y = c(ycoords1)
	xy_raster$cell_ids = c(cell_ids)
	xy_raster = adf(xy_raster)
	coordinates(xy_raster) = ~x+y
	
	
	grd_final = SpatialPixelsDataFrame(points=coordinates(xy_raster), data=adf(xy_raster$cell_ids), proj4string=grdproj)

	return(grd_final)
	}










#######################################################
# make_weeks_list:
#######################################################
#' Make a list of numbered weeks
#'
#' Make a list of weeks for a given number of days. 
#' 
#' E.g., for 15 days,
#'
#' \code{make_weeks_list(numweeks=17)}
#'
#' ...would return:
#' 
#' week1\cr
#' week1\cr
#' week1\cr
#' week1\cr
#' week1\cr
#' week1\cr
#' week1\cr
#' week2\cr
#' week2\cr
#' week2\cr
#' week2\cr
#' week2\cr
#' week2\cr
#' week2\cr
#' week3\cr
#' week3\cr
#' week3\cr
#'
#' This is useful for categorizing e.g. daily satellite data by week.
#'
#' @param numweeks Number of weeks you would like to get the daily 
#' week categories for.  Default is 53 weeks, i.e. 53*7 = 371 days.
#' @return \code{weeks_list}, a list of the week category that each day falls into
#' @export
#'   @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite NASA2001
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # Make the weeks_list
#' weeks_list = make_weeks_list(numweeks=53)
#' weeks_list
#' 
#' # Which week category you are in
#' day=1
#' weeks_list[day]
#' day=8
#' weeks_list[day]
#' day=100
#' weeks_list[day]
#' day=366
#' weeks_list[day]
#'
make_weeks_list <- function(numweeks=53)
	{
	# Week category
	
	# Take the day
	# Count the days
	
	weeks_list = NULL
	for (i in 1:numweeks)
		{
		# Write the string for the week category for 7 days
		tmpstr = paste("week", i, sep="")
		
		# Replicate 7 times
		weekdays = rep(tmpstr, 7)
		
		# Add to the weeks list
		weeks_list = c(weeks_list, weekdays)
		}

	#weeks_list
	length(weeks_list)
	
	return(weeks_list)
	}




#######################################################
# modfns_to_ftp_download_cmd: 
#######################################################
#' Make download commands for MODIS files
#'
#' Take some MODIS product filenames and produce the appropriate FTP "\code{get}" command.
#'
#' E.g., an input of:
#'
#' \code{MYD03.A2007019.1935.005.2009277181756.hdf}
#' 
#' ...is converted to:
#'
#' \code{get allData/5/MYD03/2007/019/MYD03.A2007019.1935.005.2009277181756.hdf MYD03.A2007019.1935.005.2009277181756.hdf}
#'
#' The user must specify \code{ftp_prefix}, e.g. for the products I have worked with,
#' it is as follows:
#' 
#' \code{MOD03:     get allData/5/}\cr
#' \code{MOD35_L2:  get allData/5/}\cr
#' \code{MYD03:     get allData/5/}\cr
#' \code{MYD35_L2:  get allData/5/}\cr
#' \code{MYD06_L2:  get allData/51/}\cr
#' \code{MYD06_L2:  get allData/51/}\cr
#' 
#' @param mod03_fns_dropped A list of MODIS product files (e.g. the files returned as mismatches by \code{\link{check_for_matching_geolocation_files}}
#' @param ftp_prefix The prefix that the user would like to pre-pend to the FTP directory
#' @param output_to If "screen", the FTP commands are printed to screen (and returned as a list). If "data.frame", the FTP commands are only returned as a list.  Any other string is interpreted as a filename, and the FTP commands will be saved to that as a textfile (and returned as a list).
#' @return \code{ftp_cmds}, A list of the FTP commands.
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @seealso \code{\link{check_for_matching_geolocation_files}}
#' 
modfns_to_ftp_download_cmd <- function(mod03_fns_dropped, ftp_prefix="get allData/5/", output_to="screen")
	{
	# ftp_prefix for each product:
	# MOD03:    	"get allData/5/"
	# MOD35_L2:		"get allData/5/"
	# MYD03:    	"get allData/5/"
	# MYD35_L2:		"get allData/5/"
	# MYD06_L2:		"get allData/51/"
	# MYD06_L2:		"get allData/51/"
	
	ftp_cmds = rep("", times=length(mod03_fns_dropped))
	
	# Go through each fn, and 
	for (i in 1:length(mod03_fns_dropped))
		{
		fn = mod03_fns_dropped[i]
		
		# Remove directory
		words = strsplit(fn, "/")[[1]]
		fn2 = words[length(words)]
		
		# Parse fn
		words = strsplit(fn2, "\\.")[[1]]
		
		prodname = words[1]
		A_year_day = words[2]
		
		year = substr(A_year_day, start=2, stop=5)
		day = substr(A_year_day, start=6, stop=8)
		
		
		ftp_cmds[i] = paste(ftp_prefix, prodname, "/", year, "/", day, "/", fn2, " ", fn2, sep="")
		
		}
	
	if (output_to == "data.frame")
		{
		return(ftp_cmds)
		}
	
	if (output_to == "screen")
		{
		cat("\n\nFTP cmds:\n\n", sep="")
		for (i in 1:length(mod03_fns_dropped))
			{
			cat(ftp_cmds[i], "\n", sep="")
			}
		return(ftp_cmds)
		} else {

		# output to named file
		write.table(x=ftp_cmds, file=output_to, append=FALSE, quote=FALSE, sep="\n", row.names=FALSE, col.names=FALSE)

		return(ftp_cmds)
		}
	return(ftp_cmds)
	}




#######################################################
# numslist_to_grd:
#######################################################
#' Convert a list of numbers to a grid
#'
#' Take a cloudcount raster and a number-of-observations raster and make a fraction cloud raster.
#' 
#' @param numslist The list of numbers, each corresponding to a pixel
#' @param grd A grid object, from which the projection and extend of the new grid will be taken
#' @param ydim_new Number of rows (number of pixels in height)
#' @param xdim_new Number of columns (number of pixels in width)
#' @return \code{grd_final}, as SpatialPixelsDataFrame object
#' @export
#' @seealso \code{\link{make_cloudcount_product}}
#' @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite Ackerman2010
#'   @cite GoldsmithMatzkeDawson2013
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#'
#' #######################################################
#' # Load some bit TIFs (TIFs with just the cloud indicators extracted)
#' # and add up the number of cloudy days, out of the total
#' # number of observation attempts
#' #######################################################
#' # Code excluded from CRAN check because it depends on modiscdata
#' \dontrun{
#' library(devtools)
#' # The modiscdata (MODIS c=cloud data=data) package is too big for CRAN (60 MB); so it is available on github:
#' # https://github.com/nmatzke/modiscdata
#' # If we can't get install_github() to work, try install_url():
#' # install_github(repo="modiscdata", username="nnmatzke")
#' install_url(url="https://github.com/nmatzke/modiscdata/archive/master.zip")
#' library(modiscdata)
#' tifsdir = system.file("extdata/2002_bit/", package="modiscdata")
#' tiffns = list.files(tifsdir, pattern=".tif", full.names=TRUE)
#' tiffns
#' 
#' library(rgdal)	# for readGDAL
#' 
#' # numpixels in subset
#' xdim = 538
#' ydim = 538
#' 
#' # Read the grid and the grid metadata		
#' coarsen_amount = 1
#' xdim_new = xdim / floor(coarsen_amount)
#' ydim_new = ydim / floor(coarsen_amount)
#' 
#' 
#' sum_nums = NULL
#' for (j in 1:length(tiffns))
#' 	{
#' 	fn = tiffns[j]
#' 	
#' 	grd = readGDAL(fn, output.dim=c(ydim_new, xdim_new))
#' 	
#' 	grdr = raster(grd)
#' 	pointscount_on_SGDF_points = coordinates(grd)
#' 	grdr_vals = extract(grdr, pointscount_on_SGDF_points)
#' 	
#' 	# Convert to 1/0 cloudy/not
#' 	data_grdr = grdr_vals
#' 	data_grdr[grdr_vals > 0] = 1
#' 	
#' 	grdr_cloudy = grdr_vals
#' 	grdr_cloudy[grdr_vals < 4] = 0
#' 	grdr_cloudy[grdr_vals == 4] = 1
#' 	
#' # Note: Don't run the double-commented lines unless you want to collapse different bit values.
#' 	# grdr_clear = grdr_vals
#' 	# grdr_clear[grdr_vals == 4] = 0
#' 	# grdr_clear[grdr_vals == 3] = 1
#' 	# grdr_clear[grdr_vals == 2] = 1
#' 	# grdr_clear[grdr_vals == 1] = 1
#' 	# grdr_clear[grdr_vals == 0] = 0
#' 	# 
#' 	
#' if (j == 1)
#' 	{
#' 	sum_cloudy = grdr_cloudy
#' 	#sum_not_cloudy = grdr_clear
#' 	sum_data = data_grdr
#' 	} else {
#' 	sum_cloudy = sum_cloudy + grdr_cloudy
#' 	sum_data = sum_data + data_grdr
#' 	}
#' 			
#' 	}
#' 
#'  
#' ## Calculate percentage cloudy
#' sum_nums = sum_cloudy / sum_data
#' 
#' grd_final = numslist_to_grd(numslist=sum_nums, grd=grd, ydim_new=ydim_new, xdim_new=xdim_new)
#' 
#' # Display the image (this is just the sum of a few images)
#' image(grd_final)
#' 
#' }
#' 
numslist_to_grd <- function(numslist, grd, ydim_new, xdim_new)
	{
	# Get the grd projection
	grdproj = CRS(proj4string(grd))

	tmpext2 = extent(grd)
	numcells_long = xdim_new
	numcells_lat = ydim_new
	
	xmin = attr(tmpext2, "xmin")
	xmax = attr(tmpext2, "xmax")
	ymin = attr(tmpext2, "ymin")
	ymax = attr(tmpext2, "ymax")
	xcoords = seq(from=xmin, to=xmax, by=((xmax-xmin)/numcells_long))
	ycoords = seq(from=ymax, to=ymin, by=(-1*(ymax-ymin)/numcells_lat))
	xcoords = xcoords[1:(length(xcoords)-1)]
	ycoords = ycoords[1:(length(ycoords)-1)]



	
	# Directions from here:
	# http://r-sig-geo.2731867.n2.nabble.com/Need-fast-method-to-find-indices-of-cells-in-a-3-D-grid-that-contain-x-y-z-coordinates-td2763738.html
	
	# Number the cells in the 3-D array.
	numcols = length(xcoords)
	numrows = length(ycoords)
	#cell_ids <- array(1:(numcols*numrows), dim = c(numrows,numcols))
	cell_ids <- array(numslist, dim = c(numrows,numcols))
	dim(cell_ids)
	#cell_xs <- seq(from=xmin, to=xmax, by=((xmax-xmin)/numcells_long))
	#cell_ys <- seq(from=ymax, to=ymin, by=((ymax-ymin)/numcells_long))
	
	# Get the numbers of the grid cells containing the coordinates x, y, z.
	# x, y, z are the coordinates I want to look up.  
	# xnb, ynb, and znb are the boundaries of the cells in the x, y, and z
	# directions, respectively.  
	#col_ind <- findInterval(xy_points$x, xcoords, all.inside = T)
	#row_ind <- findInterval(xy_points$y, ycoords, all.inside = T)
	#index.array <- array(c(row_ind, col_ind), dim=c(length(xy_points$x), 2))
	#indices_you_want = cell_ids[index.array]  # the numbers of the cells I want 
	
	
	xcoords1 = matrix(data=xcoords, nrow=numrows, ncol=numcols, byrow=FALSE)
	ycoords1 = matrix(data=ycoords, nrow=numrows, ncol=numcols, byrow=TRUE)
	#xcoords1
	#ycoords1
	
	xy_raster = NULL
	xy_raster$x = c(xcoords1)
	xy_raster$y = c(ycoords1)
	xy_raster$cell_ids = c(cell_ids)
	xy_raster = adf(xy_raster)
	coordinates(xy_raster) = ~x+y
	
	
	grd_final = SpatialPixelsDataFrame(points=coordinates(xy_raster), data=adf(xy_raster$cell_ids), proj4string=grdproj)

	return(grd_final)	
	}








#######################################################
# sum_bitgrid: 
#######################################################
#' Take a series of byte tifs, extract the values of a bit, and add them up
#'
#' This function is useful if you have a series of cloud product byte tifs, and need to extract 
#' a bit from each, and then add up (e.g. to get the total number of cloudy observations).
#' 
#' @param fns A list of filenames which can be auto-read by readGDAL (e.g. geoTIFFs)
#' @param numits Number of iterations to run for (default = length(fns))
#' @param bitnums The bit(s) you would like to extract. If length(bitnums) = 1, do default behavior.  If the length is 2, then run the alternative function.
#' @param ydim_new The dimensions of the extracted image.  This is REQUIRED, since each day of a MODIS image may have different pixel dimensions over the same location. readGDAL will sample each image (nearest neighbor pixel) to the correct dimensions.
#' @param xdim_new The dimensions of the extracted image.  This is REQUIRED, since each day of a MODIS image may have different pixel dimensions over the same location. readGDAL will sample each image (nearest neighbor pixel) to the correct dimensions.
#' @param test_for_1 What your extracted bit(s) should match. If 1 (default), just count up 1s.  If a string, e.g. "!=11", use that in an if-then statement to get the 1s to count
#' @return \code{grdr_vals_bits}, a matrix of 0/1 values for the bit in question
#' @export
#' @seealso \code{\link{extract_bit}}
#' @seealso \code{\link{numslist_to_grd}}
#' @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite Ackerman2010
#'   @cite GoldsmithMatzkeDawson2013
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' 
#' # See numslist_to_grd example for an idea of how this works.
#'
sum_bitgrid <- function(fns, numits=length(fns), bitnums, ydim_new, xdim_new, test_for_1=1)
	{
	require(rgdal)	# for readGDAL
	
	# Numits is a useful shortcut if you want to run fewer files, but don't want
	# subset each list of files; however, if 
	# Error catch
	if (numits > length(fns))
		{
		numits = length(fns)
		}
	
	for (i in 1:numits)
		{
		############################
		# Look for bit 0
		############################
		fn = fns[i]
		#grd = readGDAL(fn)
		grd = rgdal::readGDAL(fn, output.dim=c(ydim_new, xdim_new))
		#print(summary(grd))
		#image(grd)
	
		# Different behavior, if this is 1 or 2 bits
		if (length(bitnums) == 1)
			{
			# Extract just bit 0	
			grdr_vals_bits = get_bitgrid(grd, bitnum=bitnums[1])
			
			# If you want something other than "1" to equal 1
			if (test_for_1 != 1)
				{
				# Make a new grid and set to 0
				grdr_vals_bits2 = grdr_vals_bits
				grdr_vals_bits2[!is.na(grdr_vals_bits2)] = 0
				
				# Put the 1s into the matching pixels
				cmdstr = paste("grdr_vals_bits2[grdr_vals_bits", test_for_1, "] = 1", sep="")
				eval(parse(text=cmdstr))
				
				grdr_vals_bits = grdr_vals_bits2
				}
			}
		# If you want a string of 2 bits in each pixel
		if (length(bitnums) == 2)
			{
			# Extract a 2-digit string representing the 2 bits	
			grdr_vals_bitstrings = get_bitgrid_2bits(grd, bitnums=bitnums)
			
			# If you want something other than "1" to equal 1, which you 
			# SHOULD in this case, since we are dealing with a string
			#
			# e.g., a common string would be "!= '11'", because of 
			# 11 = confident clear, anything else is possible or probable cloud
			#
			if (test_for_1 == 1)
				{
				print("ERROR: Stop -- when you have a 2-bit string per pixel, you need a 2-bit search term to count something, e.g. != '11'")
				}
			
			if (test_for_1 != 1)
				{
				# Make a new grid and set to 0
				grdr_vals_bits2 = grdr_vals_bitstrings
				grdr_vals_bits2[!is.na(grdr_vals_bits2)] = as.numeric(0)
				grdr_vals_bits2 = as.numeric(grdr_vals_bits2)
				
				# Put the 1s into the matching pixels
				cmdstr = paste("grdr_vals_bits2[grdr_vals_bitstrings", test_for_1, "] = 1", sep="")
				eval(parse(text=cmdstr))
				
				grdr_vals_bits = grdr_vals_bits2
				}
			
			}
			
		# Add up the values...
		# If i=1, set up a master image to do the adding
		if (i == 1)
			{
			# Set up the gridpixels list, set all values to 0
			master_bitgridvals = grdr_vals_bits
			master_bitgridvals[master_bitgridvals > -100] = 0
			master_bitgridvals = master_bitgridvals + grdr_vals_bits
			}
		else
			{
			# Add 'em up!	
			master_bitgridvals = master_bitgridvals + grdr_vals_bits		
			}
	
			
		}

	print("Unique grid values:")
	print(sort(unique(master_bitgridvals)))

	return(master_bitgridvals)
	}








#######################################################
# write_MRTSwath_param_file:
#######################################################
#' Write a parameter control file for MRTSwath
#'
#' MRTSwath is the "MODIS Reprojection Tool for swath products".  See:
#' \url{https://lpdaac.usgs.gov/tools/modis_reprojection_tool_swath}).
#' 
#' If you want this function to use MRTSwath tool successfully, you should 
#' add the directory with the MRTSwath executable to the default R PATH
#' by editing \code{~/.Rprofile}.
#'
#' This function hard-codes these options into the parameter file:\cr
#' * all the bands are extracted\cr
#' * the output file is a GeoTIFF\cr
#' * the output projection is Geographic (plain unprojected Latitude/Longitude)\cr
#' * the resampling is Nearest Neighbor (NN), which of course is the only one which makes sense when the pixels encode bytes that encode bits that encode discrete classification results, 0/1 error flags, etc.\cr
#'
#' MRTswath can do many other projections and output formats; users can modify this function to run those options.
#'
#' @param prmfn The name of the parameter/control file which will be the input to MRTSwath's \code{swath2grid} function.
#' @param tifsdir The directory to save the output TIF files in
#' @param modfn The filename of the MODIS data
#' @param geoloc_fn The filename of the corresponding geolocation file (annoyingly, this is a much larger
#' file than the data file!)
#' @param ul_lon Upper left (ul) longitude (x-coordinate) for subsetting
#' @param ul_lat Upper left (ul) latitude (y-coordinate) for subsetting
#' @param lr_lon Lower right (lr) longitude (x-coordinate) for subsetting
#' @param lr_lat Lower right (lr) latitude (y-coordinate) for subsetting
#' @return \code{prmfn} The name of the temporary parameter file
#' @export
#' @seealso \code{\link{run_swath2grid}}
#' @seealso \url{http://landweb.nascom.nasa.gov/cgi-bin/QA_WWW/newPage.cgi?fileName=hdf_filename}
#'   @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite NASA2001
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#'
#' # Source MODIS files (both data and geolocation)
#' # Code excluded from CRAN check because it depends on modiscdata
#' \dontrun{
#' library(devtools)
#' # The modiscdata (MODIS c=cloud data=data) package is too big for CRAN (60 MB); so it is available on github:
#' # https://github.com/nmatzke/modiscdata
#' # If we can't get install_github() to work, try install_url():
#' # install_github(repo="modiscdata", username="nnmatzke")
#' install_url(url="https://github.com/nmatzke/modiscdata/archive/master.zip")
#' library(modiscdata)
#' moddir = system.file("extdata/2002raw/", package="modiscdata")
#' 
#' # Get the matching data/geolocation file pairs
#' fns_df = check_for_matching_geolocation_files(moddir, modtxt="MYD35_L2", geoloctxt="MYD03")
#' fns_df
#' 
#' # Resulting TIF files go in this directory
#' tifsdir = getwd()
#' 
#' 
#' # Box to subset
#' ul_lat = 13
#' ul_lon = -87
#' lr_lat = 8
#' lr_lon = -82
#' 
#' for (i in 1:nrow(fns_df))
#' 	{
#' 	
#' 	prmfn = write_MRTSwath_param_file(prmfn="tmpMRTparams.prm", tifsdir=tifsdir, modfn=fns_df$mod35_L2_fns[i], geoloc_fn=fns_df$mod03_fns[i], ul_lon=ul_lon, ul_lat=ul_lat, lr_lon=lr_lon, lr_lat=lr_lat)
#' 	print(scan(file=prmfn, what="character", sep="\n"))
#' 	
#' 	}
#' }
#' 
write_MRTSwath_param_file <- function(prmfn="tmpMRTparams.prm", tifsdir, modfn, geoloc_fn, ul_lon, ul_lat, lr_lon, lr_lat)
	{
	# Initialize the list of lines in the parameter file
	prmfile = NULL
	pnum = 0
	prmfile[[(pnum=pnum+1)]] = 	" "

	# Input files
	prmfile[[(pnum=pnum+1)]] = 	paste("INPUT_FILENAME = ", modfn, sep="")
	prmfile[[(pnum=pnum+1)]] = 	" "
	prmfile[[(pnum=pnum+1)]] = 	paste("GEOLOCATION_FILENAME = ", geoloc_fn, sep="")
	prmfile[[(pnum=pnum+1)]] = 	" "
	prmfile[[(pnum=pnum+1)]] = 	paste("INPUT_SDS_NAME = Cloud_Mask, 1, 1, 1, 1, 1, 1", sep="")
	
	# Subset parameters
	prmfile[[(pnum=pnum+1)]] = 	" "
	prmfile[[(pnum=pnum+1)]] = 	paste("OUTPUT_SPATIAL_SUBSET_TYPE = LAT_LONG", sep="")
	prmfile[[(pnum=pnum+1)]] = 	paste("OUTPUT_SPACE_UPPER_LEFT_CORNER (LONG LAT) =", ul_lon, ul_lat, sep=" ")
	prmfile[[(pnum=pnum+1)]] = 	paste("OUTPUT_SPACE_LOWER_RIGHT_CORNER (LONG LAT) = ", lr_lon, lr_lat, sep=" ")
	
	
	# Output filename
	prmfile[[(pnum=pnum+1)]] = 	" "
	
	outfn = gsub(pattern=".hdf", replacement=".tif", modfn)
	outfn = extract_fn_from_path(fn_with_path=outfn)
	outfn = slashslash(paste(tifsdir, outfn, sep="/"))

	prmfile[[(pnum=pnum+1)]] = 	paste("OUTPUT_FILENAME = ", outfn, sep="")
	prmfile[[(pnum=pnum+1)]] = 	paste("OUTPUT_FILE_FORMAT = GEOTIFF_FMT", sep="")

	# Reprojection information (for Geographic Projection, with nearest-neighbor resampling)
	prmfile[[(pnum=pnum+1)]] = 	" "
	prmfile[[(pnum=pnum+1)]] = 	"KERNEL_TYPE (CC/BI/NN) = NN"
	prmfile[[(pnum=pnum+1)]] = 	" "
	prmfile[[(pnum=pnum+1)]] = 	"OUTPUT_PROJECTION_NUMBER = GEO"
	prmfile[[(pnum=pnum+1)]] = 	" "
	prmfile[[(pnum=pnum+1)]] = 	"OUTPUT_PROJECTION_PARAMETER = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"
	prmfile[[(pnum=pnum+1)]] = 	" "
	prmfile[[(pnum=pnum+1)]] = 	"OUTPUT_PROJECTION_SPHERE = 8"
	prmfile[[(pnum=pnum+1)]] = 	" "

	#prmfile
	
	
	# Write the MRTSwath tool parameter file
	write.table(x=prmfile, file=prmfn, append=FALSE, quote=FALSE, sep="\n", row.names=FALSE, col.names=FALSE)
	#moref(prmfn)
		
	return(prmfn)
	}




#######################################################
# run_swath2grid:
#######################################################
#' Run MRTSwath swath2grid tool
#'
#' MRTSwath is the "MODIS Reprojection Tool for swath products".  See:
#' \url{https://lpdaac.usgs.gov/tools/modis_reprojection_tool_swath}).
#' 
#' If you want this function to use MRTSwath tool successfully, you should 
#' add the directory with the MRTSwath executable to the default R PATH
#' by editing \code{~/.Rprofile}.
#'
#' @param mrtpath This is the path to the MRTSwath executable \code{swath2grid}. If your \code{~/.Rprofile}
#' file has the location of \code{swath2grid} in the PATH, then you can just use \code{mrtpath="swath2grid"}.
#' Otherwise, the user must provide the full path to swath2grid.
#' @param prmfn The name of the parameter/control file which will be the input to MRTSwath's \code{swath2grid} function.
#' @param tifsdir The directory to save the output TIF files in
#' @param modfn The filename of the MODIS data
#' @param geoloc_fn The filename of the corresponding geolocation file (annoyingly, this is a much larger
#' file than the data file!)
#' @param ul_lon Upper left (ul) longitude (x-coordinate) for subsetting
#' @param ul_lat Upper left (ul) latitude (y-coordinate) for subsetting
#' @param lr_lon Lower right (lr) longitude (x-coordinate) for subsetting
#' @param lr_lat Lower right (lr) latitude (y-coordinate) for subsetting
#' @return \code{cmdstr} The string giving the system command that ran \code{swath2grid}
#' @export
#' @seealso \code{\link{write_MRTSwath_param_file}}
#' @seealso \url{http://landweb.nascom.nasa.gov/cgi-bin/QA_WWW/newPage.cgi?fileName=hdf_filename}
#'   @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite NASA2001
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' #######################################################
#' # Run MRTSwath tool "swath2grid"
#' #######################################################
#' 
#' # Source MODIS files (both data and geolocation)
#' # Code excluded from CRAN check because it depends on modiscdata
#' \dontrun{
#' library(devtools)
#' # The modiscdata (MODIS c=cloud data=data) package is too big for CRAN (60 MB); so it is available on github:
#' # https://github.com/nmatzke/modiscdata
#' # If we can't get install_github() to work, try install_url():
#' # install_github(repo="modiscdata", username="nnmatzke")
#' install_url(url="https://github.com/nmatzke/modiscdata/archive/master.zip")
#' library(modiscdata)
#' moddir = system.file("extdata/2002raw/", package="modiscdata")
#' 
#' # Get the matching data/geolocation file pairs
#' fns_df = check_for_matching_geolocation_files(moddir, modtxt="MYD35_L2", geoloctxt="MYD03")
#' fns_df
#' 
#' # Resulting TIF files go in this directory
#' tifsdir = getwd()
#' 
#' 
#' # Box to subset
#' ul_lat = 13
#' ul_lon = -87
#' lr_lat = 8
#' lr_lon = -82
#' 
#' for (i in 1:nrow(fns_df))
#' 	{
#' 	
#'	prmfn = write_MRTSwath_param_file(prmfn="tmpMRTparams.prm", tifsdir=tifsdir, modfn=fns_df$mod35_L2_fns[i], geoloc_fn=fns_df$mod03_fns[i], ul_lon=ul_lon, ul_lat=ul_lat, lr_lon=lr_lon, lr_lat=lr_lat)
#'	print(scan(file=prmfn, what="character", sep="\n"))
#' 	
#'	run_swath2grid(mrtpath="swath2grid", prmfn="tmpMRTparams.prm", tifsdir=tifsdir, modfn=fns_df$mod35_L2_fns[i], geoloc_fn=fns_df$mod03_fns[i], ul_lon=ul_lon, ul_lat=ul_lat, lr_lon=lr_lon, lr_lat=lr_lat)
#'	
#' 	}
#'
#' list.files(tifsdir, pattern=".tif", full.names=TRUE)
#' }
#'
run_swath2grid <- function(mrtpath="swath2grid", prmfn="tmpMRTparams.prm", tifsdir, modfn, geoloc_fn, ul_lon, ul_lat, lr_lon, lr_lat)
	{
	# Check for the existence of mrtpath
	which_result = system(paste("which ", mrtpath, sep=""), intern=TRUE)
	which_result
	
	if (length(which_result) == 0)
		{
		return(paste("Error: mrtpath (", mrtpath, ") does not correspond to a file. Having swath2grid from MRTswatch is required for this function.", sep=""))
		}
	
	# Write the temporary parameter file
	prmfn = write_MRTSwath_param_file(prmfn=prmfn, tifsdir=tifsdir, modfn=modfn, geoloc_fn=geoloc_fn, ul_lon=ul_lon, ul_lat=ul_lat, lr_lon=lr_lon, lr_lat=lr_lat)
	
	# Run MRTSwath tool (swath2grid)
	cmdstr = paste(mrtpath, " -pf=", prmfn, sep="")
	system(cmdstr)
	
	return(cmdstr)
	}





#######################################################
# yearday_to_date:
#######################################################
#' Convert a year + a day number to a date
#'
#' The filenames of MODIS images contain the following date
#' information: MOD35_L2.Ayyyyddd.hhhh.etc.
#'
#'     MODLAND Level 2 products
#'     ESDT.AYYYYDDD.HHMM.CCC.YYYYDDDHHMMSS.hdf
#'     
#'     ESDT = Earth Science Data Type name (e.g., MOD14)
#'     YYYYDDD = MODIS acquisition year and Julian day
#'     HHMM = MODIS acquisition UTC time
#'     CCC = Collection number
#'     YYYYDDDHHMMSS = Processing Year, Julian day and UTC Time
#'     hdf = Suffix denoting HDF file 
#'
#' DDD is the day of the year, from 001 to 365 (or 366 for leap years)
#' 
#' This is mildly annoying to interpret as e.g. months for graphing
#' cloudy days per month, so \code{\link{yearday_to_date}} converts a (numeric)
#' year and day to the calendar date.
#' 
#' @param year The year, read as numeric
#' @param day The day of the year, read as numeric, from 1 to 366
#' @return \code{newdate} A list with three items: \code{$month}, \code{$day}, \code{$year}
#' @export
#' @seealso \code{\link{make_POSIXct_date}}
#' @seealso \code{\link{yearday_to_date}}
#' @seealso \url{http://landweb.nascom.nasa.gov/cgi-bin/QA_WWW/newPage.cgi?fileName=hdf_filename}
#'   @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite NASA2001
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' yearday_to_date(year=2012, day=364)
#' # $month
#' # [1] 12
#' # 
#' # $day
#' # [1] 29
#' # 
#' # $year
#' # [1] 2008
#'
yearday_to_date <- function(year=2012, day=1)
	{
	# The "date" package is required
	# install.packages("date")
	require(date)
	
	# Force conversion to numeric
	year = as.numeric(year)
	day = as.numeric(day)
	
	# the date package calculates Julian days from January 1, 1960
	# Start by getting the Julian day of Jan. 1 of the year of interest
	days_since_1960_for_Jan1 = as.numeric(mdy.date(month=1, day=1, year=year))
	days_since_1960_for_Jan1
	
	# Add the day number (minus 1 for January 1)
	days_since_1960_for_image = days_since_1960_for_Jan1 - 1 + day
	
	# Convert to 
	newdate = date.mdy(sdate=days_since_1960_for_image)
	
	return(newdate)
	}


#######################################################
# yearmonthday_to_julianday:
#######################################################
#' Get the julian day for a year/month/day date
#'
#' This function uses the date package's \code{\link[date]{mdy.date}} function to get the 
#' Julian day (count of the day in the year, from 1-365, or 1-366 for leap years) from
#' the input year, month, and day.
#' 
#' @param year as numeric, 4 digits
#' @param month as numeric, 1-2 digits
#' @param day as numeric, 1-2 digits
#' @return julian_day, the julian day between 1-365 or 1-366 (for leap years)
#' @export
#' @seealso \code{\link{yearday_to_date}}
#' @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite Ackerman2010
#'   @cite GoldsmithMatzkeDawson2013
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' year=2011; month=06; day=15
#' yearmonthday_to_julianday(year, month, day)
#'
yearmonthday_to_julianday <- function(year, month, day)
	{
	require(date)
	
	# Get the Julian day for this day this year
	days_since_1960_for_Jan1 = as.numeric(mdy.date(month=1, day=1, year=year))
	days_since_1960_for_current_date = as.numeric(mdy.date(month=month, day=day, year=year))

	# Add the day number (plus 1 for January 1)
	julian_day = days_since_1960_for_current_date - days_since_1960_for_Jan1 + 1

	return(julian_day)
	}






