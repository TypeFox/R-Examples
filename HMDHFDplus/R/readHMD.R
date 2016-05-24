
# FILE CONTENTS:

# 1) readHMD()
# 2) readHMDweb()
# 3) readJMDweb()

############################################################################
# readHMD()
############################################################################

#'
#' @title \code{readHMD()} reads a standard HMD .txt table as a \code{data.frame}
#' 
#' @description This calls \code{read.table()} with all the necessary defaults to avoid annoying surprises. The Age column is also stripped of \code{"+"} and converted to integer, and a logical indicator column called \code{OpenInterval} is added to show where these were located. If the file contains population counts, values are split into two columns for Jan 1 and Dec 31 of the year. Output is invisibly returned, so you must assign it to take a look. This is to avoid lengthy console printouts. 
#' 
#' @param filepath path or connection to the HMD text file, including .txt suffix.
#' @param ... other arguments passed to \code{read.table}, not likely needed.
#' @param fixup logical. Should columns be made more user-friendly, e.g., forcing Age to be integer?
#' 
#' @return data.frame of standard HMD output, except the Age column has been cleaned, and a new open age indicator column has been added. If the file is Population.txt or Population5.txt, there will be two columns each for males and females.
#' 
#' @details Population counts in the HMD typically refer to Jan 1st. One exception are years in which a territorial adjustment has been accounted for in estimates. For such years, `YYYY-` refers to Dec 31 of the year before the adjustment, and `YYYY+` refers to Jan 1 directly after the adjustment (adjustments are always made Jan 1st). In the data, it will just look like two different estimates for the same year, but in fact it is a definitional change or similar. In order to remove headaches from potential territorial adjustments in the data, we simply create two columns, one for January 1st (e.g.,\code{"Female1"}) and another for Dec 31st (e.g.,\code{"Female2"}) . One can recover the adjustment coefficient for each year by taking the ratio $$Vx = P1(t+1) / P2(t)$$. In most years this will be 1, but in adjustment years there is a difference. This must always be accounted for when calculating rates and exposures. Argument \code{fixup} is outsourced to \code{HMDparse()}.
#' 
#' @importFrom utils read.table
#' 
#' @export
#' 
#' @note function written by Tim Riffe.
#' 
readHMD <- function(filepath, fixup = TRUE, ...){
  DF              <- read.table(file = filepath, header = TRUE, skip = 2, na.strings = ".", as.is = TRUE, ...)
  if (fixup){
    DF        <- HMDparse(DF, filepath)
  }
  invisible(DF)
}

############################################################################
# readHMDweb()
############################################################################

#'
#' @title readHMDweb a basic HMD data grabber.
#' 
#' @description This is a basic HMD data grabber, based on Carl Boe's original \code{HMD2R()}. It will only grab a single HMD statistical product from a single country. Some typical R pitfalls are removed: The Age column is coerced to integer, while an AgeInterval column is created. Also Population counts are placed into two columns, for Jan. 1st and Dec. 31 of the same year, so as to remove headaches from population universe adjustments, such as territorial changes. Fewer options means less to break. To do more sophisticated data extraction, iterate over country codes or statistical items. Reformatting can be done outside this function using, e.g., \code{long2mat()}. Argument \code{fixup} is outsourced to \code{HMDparse()}.
#'
#' @param CNTRY HMD population letter code. If not spelled right, or not specified, the function provides a selection list. Only 1.
#' @param item the statistical product you want, e.g., \code{"fltper_1x1"}. Only 1.
#' @param username usually the email address you registered with the HMD under. If left blank, you'll be prompted. Do that if you don't mind the typing and prefer not to save your username in your code.
#' @param password Your HMD password. If left blank, you'll be prompted. Do that if you don't mind the typing and prefer not to save your password in your code.
#' @param fixup logical. Should columns be made more user-friendly, e.g., forcing Age to be integer?
#' 
#' @return data.frame of the HMD product, read as as \code{readHMD()} would read it.
#'
#' @details You need to register for HMD: \url{www.mortality.org}. It is advised to pass in your credentials as named vectors rather than directly as character strings, so that they are not saved directly in your code. See examples. One option is to just save them in your Rprofile file.
#' 
#' @importFrom RCurl getURL
#' @importFrom RCurl getCurlHandle
#' @importFrom RCurl getCurlInfo
#' @importFrom RCurl url.exists
#' @importFrom utils read.csv
#' @importFrom utils read.table
#' 
#' @export
#' 
readHMDweb <- function(CNTRY = NULL, item = NULL, username = NULL, password = NULL, fixup = TRUE){
	## based on Carl Boe's RCurl tips
	# modified by Tim Riffe 
	
	# let user input name and password
	if (is.null(username)){
		if (interactive()){
			cat("\ntype in HMD username (usually your email, quotes not necessary):\n")
			username <- userInput(FALSE)
		} else {
			stop("if username and password not given as arguments, the R session must be interactive.")
		}
	}
	if (is.null(password)){
		if (interactive()){
			cat("\ntype in HMD password:\n")
			password <-  userInput(FALSE)
		} else {
			stop("if username and password not given as arguments, the R session must be interactive.")
		}
	}
	
	urlbase         <- "http://www.mortality.org/hmd"
#    tf <- tempfile()
#    on.exit(unlink(tf))
	this.url    <- "http://www.mortality.org/countries.csv"
	cntries     <- RCurl::getURL(this.url)
	ctrylist    <- read.csv(text = cntries,header=TRUE,as.is=TRUE);
	ctrylookup  <- data.frame(Country=ctrylist$Country, CNTRY=ctrylist$Subpop.Code.1, stringsAsFactors = FALSE)
	
	# get CNTRY
	if(is.null(CNTRY)){    
		cat("\nCNTRY missing\n")
		if (interactive()){
			CNTRY <- select.list(choices = ctrylookup$CNTRY, multiple = FALSE, title = "Select Country Code")
		} else {
			stop("CNTRY should be one of these:\n",paste(ctrylookup$CNTRY, collapse = ",\n"))
		}
	}
	if (!(CNTRY %in% ctrylookup$CNTRY)){
		cat("\nCNTRY not found\n")
		if (interactive()){
			CNTRY <- select.list(choices = ctrylookup$CNTRY, multiple = FALSE, title = "Select Country Code")
		} else {
			stop("CNTRY should be one of these:\n",paste(ctrylookup$CNTRY, collapse = ",\n"))
		}
	}
	stopifnot(length(CNTRY) == 1)
	
	this.pw <- paste(username, password, sep = ":")
	
	## reuse handle, reduce connection starts
	handle <- RCurl::getCurlHandle(userpwd = this.pw)
	
	dirjunk <- RCurl::getURL(file.path("www.mortality.org", "hmd", CNTRY,
					paste0("STATS",.Platform$file.sep)), curl = handle)
	
	if (RCurl::getCurlInfo(handle)$response.code == 401) {
		stop("Authentication rejected: please check your username and password")
	}
	dirjunk <- RCurl::getURL(file.path("www.mortality.org","hmd",CNTRY,"STATS/"), curl=handle)
	
	# check if authentication fails
	if (RCurl::getCurlInfo(handle)$response.code == 401){
		stop("Authentication rejected: please check your username and password")
	}
	# sometime redirects will break this, so we do it manually if necessary...
	if (RCurl::getCurlInfo(handle)$response.code == 301){
		dirjunk <- RCurl::getURL(getCurlInfo(handle)$redirect.url, curl = handle)
	}
	
	# TR: this is the kind of parsing I hate. Gotta be a better way out there.
	parts <- gsub(pattern = "\\\"",
			replacement = "",
			unlist(lapply(strsplit(unlist(strsplit(dirjunk
													,split="href=")),
									split = ">"),"[[",1)))
	allitems <- gsub(pattern = ".txt",replacement = "",parts[grepl(parts,pattern=".txt")])
	
	if (is.null(item)){
		if (interactive()){
			cat("\nThe following items are available for", CNTRY,"\n")
			item <- select.list(choices = allitems, 
					multiple = FALSE,
					title = "Select one")
		} else {
			stop("item must be one of the following for",CNTRY,paste(allitems,collapse=",\n"))
		}
	}
	if (!item %in% allitems){
		if (interactive()){
			if (any(grepl(allitems, pattern = item))){
				cat("\nMust specify item fully\n")    
				item <- select.list(choices = allitems[grepl(allitems, pattern = item)], 
						multiple = FALSE,
						title = "Select one")
			} else {
				cat("\nThe following items are available for", CNTRY,"\n")
				item <- select.list(choices = allitems, 
						multiple = FALSE,
						title = "Select one")
			}
		} else {
			stop("item must be one of the following for",CNTRY,paste(allitems,collapse=",\n"))
		}
	}
	
	# build url: 
    # TR: presumably these links are composed with the same separators everywhere?
	HMDurl <- paste("www.mortality.org", "hmd", CNTRY, "STATS", paste0(item, ".txt"), sep = "/")
	
	#check it exists:
    # TR: this is like way extra, since both CNTRY and item have gone through filters by now
	if (RCurl::url.exists(HMDurl,curl=handle)){
		handle <- RCurl::getCurlHandle(userpwd = this.pw)
		# grab the data
		dataIN  <- RCurl::getURL(HMDurl, curl=handle)
		
		# rest of this lifted from readHMD()
		DF      <- read.table(text = dataIN, header = TRUE, skip = 2, na.strings = ".", as.is = TRUE)
		if (fixup){
			DF        <- HMDparse(DF, filepath = item)
		}
		
		return(invisible(DF))
	} else {
		cat("\nSorry, something was wrong with the query\nPossibly a typo?\n")
	}
} # end readHMDweb()

############################################################################
# readJMDweb()
############################################################################

#'
#' @title read data from the Japan Mortality Database into R
#' 
#' @description JMD data are formatted exactly as HMD data. This function simply parses the necessary url together given a prefecture code and data item (same nomenclature as HMD). Data is parsed using \code{HMDparse()}, which converts columns into useful and intuitive classes, for ready-use. See \code{?HMDparse} for more information on type conversions. No authentification is required for this database. Only a single item/prefecture is downloaded. Loop for more complex calls (See examples). The prefID is not appended as a column, so be mindful of this if appending several items together into a single \code{data.frame}. Note that at the time of this writing, the finest Lexis resolution for prefectural lifetables is 5x5 (5-year, 5-year age groups). Raw data are, however, provided in 1x1 format, and deaths are also available in triangles.
#' 
#' @param prefID a single prefID 2-digit character string, ranging from \code{"00"} to \code{"47"}.
#' @param item the statistical product you want, e.g., \code{"fltper_5x5"}. Only 1.
#' @param fixup logical. Should columns be made more user-friendly, e.g., forcing Age to be integer?
#' @param ... extra arguments ultimately passed to \code{read.table()}. Not likely needed.
#' 
#' @return \code{data.frame} of the data item is invisibly returned
#' 
#' @details No details of note. This database in independently maintained, so file types/locations are subject to change. If this happens, please notify the package maintainer.
#' 
#' @importFrom RCurl url.exists
#' 
#' @export 
#' 
#' @examples 
#' \dontrun{
#' library(HMDHFDplus)
#' # grab prefecture codes (including All Japan)
#' prefectures <- getJMDprefectures()
#' # grab all mltper_5x5
#' # and stick into long data.frame: 
#' mltper <- do.call(rbind, lapply(prefectures, function(prefID){
#'                    Dat        <- readJMDweb(prefID = prefID, item = "mltper_5x5", fixup = TRUE)
#'                    Dat$PrefID <- prefID
#'                    Dat
#' }))
#' }
#' 
readJMDweb <- function(prefID = "01", item = "Deaths_5x5", fixup = TRUE, ...){
	JMDurl      <- paste("http://www.ipss.go.jp/p-toukei/JMD",
			         prefID, "STATS", paste0(item, ".txt"), sep = "/")
	if (RCurl::url.exists(JMDurl)){
		con         <- url(JMDurl)
		Dat         <- readHMD(con, fixup = fixup, ...)
		#close(con)
		return(invisible(Dat))
	} else {
		cat("Either the prefecture code or data item are not available\nCheck names.\nNULL returned\n")
		NULL
	}
}
# item <- "mltper_5x5";Dat <- readHMD(con, fixup = TRUE)
############################################################################
# readCHMDweb()
############################################################################

#'
#' @title read data from the Canadian Human Mortality Database into R
#' 
#' @description CHMD data are formatted exactly as HMD data. This function simply parses the necessary url together given a province code and data item (same nomenclature as HMD). Data is parsed using \code{HMDparse()}, which converts columns into useful and intuitive classes, for ready-use. See \code{?HMDparse} for more information on type conversions. No authentification is required for this database. Only a single item/prefecture is downloaded. Loop for more complex calls (See examples). The provID is not appended as a column, so be mindful of this if appending several items together into a single \code{data.frame}. Note that at the time of this writing, the finest Lexis resolution for prefectural lifetables is 5x5 (5-year, 5-year age groups). Raw data are, however, provided in 1x1 format, and deaths are also available in triangles. Note that cohort data are not produced for Canada at this time (but you could produce such data by starting with the \code{Deaths\_Lexis} file...).
#' 
#' @param provID a single provID 3 character string, as returned by \code{getCHMDprovinces()}.
#' @param item the statistical product you want, e.g., \code{"fltper_5x5"}. Only 1.
#' @param fixup logical. Should columns be made more user-friendly, e.g., forcing Age to be integer?
#' @param ... extra arguments ultimately passed to \code{read.table()}. Not likely needed.
#' 
#' @return \code{data.frame} of the data item is invisibly returned
#' 
#' @details This database is curated independently from the HMD/HFD family, and so file types and locations may be subject to change. If this happens, please notify the package maintainer.
#' 
#' @export 
#' 
#' @importFrom RCurl url.exists
#' 
#' @examples 
#' \dontrun{
#' library(HMDHFDplus)
#' # grab province codes (including All Canada)
#' provs <- getCHMDprovinces()
#' # grab all mltper_5x5  
#' # and stick into long data.frame: 
#' mltper <- do.call(rbind, lapply(provs, function(provID){
#'                    Dat        <- readCHMDweb(provID = provID, item = "mltper_5x5", fixup = TRUE)
#'                    Dat$provID <- provID
#'                    Dat
#' }))
#' }
#' 

readCHMDweb <- function(provID = "can", item = "Deaths_1x1", fixup = TRUE, ...){
	CHMDurl         <- paste("http://www.prdh.umontreal.ca/BDLC/data/",
			             provID, paste0(item, ".txt"), sep = "/")

	if (RCurl::url.exists(CHMDurl)){
		con         <- url(CHMDurl)
		Dat         <- readHMD(con, fixup = fixup, ...)
		#close(con)
		return(invisible(Dat))
	} else {
		cat("Either the prefecture code or data item are not available\nCheck names.\nNULL returned\n")
		NULL
	}
	
	
	invisible(Dat)
}

