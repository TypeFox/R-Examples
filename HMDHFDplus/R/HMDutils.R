
# FILE CONTENTS:

# 1) HMDparse()
# 2) getHMDcountries()
# 3) getJMDprefectures()
# 4) getCHMDprovinces()

############################################################################
# 1) HMDparse()
############################################################################
# note to self and others: if the HFD made metadata openly available, such as the
# below-used csv file, then everything would be so much easier!
############################################################################
#'
#' @title internal function for modifying freshly read HMD data in its standard form
#' 
#' @description called by \code{readHMD()} and \code{readHMDweb()}. We assume there are no factors in the given data.frame and that it has been read in from the raw text files using something like: \code{ read.table(file = filepath, header = TRUE, skip = 2, na.strings = ".", as.is = TRUE)}. This function is visible to users, but is not likely needed directly. 
#' 
#' @param DF a data.frame of HMD data, freshly read in.
#' @param filepath just to check if these are population counts from the name. 
#' 
#' @return DF same data.frame, modified so that columns are of a useful class. If there were open age categories, such as \code{"-"} or \code{"+"}, this information is stored in a new dummy column called \code{OpenInterval}.
#' 
#' @details This parse routine is based on the subjective opinions of the author...
#' 
#' @export
#'

HMDparse <- function(DF, filepath){
	if (any(grepl("age", tolower(colnames(DF))))){
		Pluses          <- grepl(pattern = "\\+", DF$Age )
		DF$Age          <- age2int(DF$Age)
		DF$OpenInterval <- Pluses
	}
	# Population.txt is a special case:
	if (grepl("pop", tolower(filepath))){
		# what years do we have?
		all.years   <- sort(unique(age2int(DF$Year)))
		# make indicators:
		Pluses      <- grepl(pattern = "\\+", DF$Year )
		Minuses     <- grepl(pattern = "\\-", DF$Year )
		# split out DF into two parts sum(Minuses) 
		Jan1i       <- DF$Year %in% as.character(all.years[-length(all.years)]) | Minuses
		Dec31i      <- DF$Year %in% as.character(all.years[-1]) | Pluses
		Jan1        <- DF[Jan1i, ]
		Dec31       <- DF[Dec31i, ]
		
		Jan1$Year   <- age2int(Jan1$Year)
		Dec31$Year  <- age2int(Dec31$Year)
		
		# now stick back together just the parts we need:
		cols1       <- match(c("female","male","total"),tolower(colnames(Jan1)))
		cols2       <- match(c("female","male","total"),tolower(colnames(Dec31)))
		colnames(Jan1)[cols1]   <- paste0(colnames(Jan1)[cols1],1)
		colnames(Dec31)[cols2]  <- paste0(colnames(Dec31)[cols2],2)
		DF          <- cbind(Jan1, Dec31[,cols2])
		# finally reorganize columns:
		orgi        <- grepl("male",tolower(colnames(DF))) | grepl("total",tolower(colnames(DF)))
		DF          <- cbind(DF[, !orgi], DF[, orgi])
	}
	
	if (any(grepl("year", tolower(colnames(DF))))){
		DF$Year          <- age2int(DF$Year)
	}
	if (any(grepl("cohort", tolower(colnames(DF))))){
		DF$Cohort          <- age2int(DF$Cohort)
	}
	invisible(DF)
}

############################################################################
# 2) getHMDcountries()
############################################################################

#' @title internal function for grabbing the HMD country short codes. 
#'
#' @description This function is called by \code{readHMDweb()} and is separated here for modularity. Assumes you have an internet connection.
#' 
#' @return a vector of HMD country short codes.
#' 
#' @importFrom utils read.csv
#' 
#' @export

getHMDcountries <- function(){
	HMDXXX  <- read.csv("http://www.mortality.org/countries.csv",stringsAsFactors = FALSE)
	HMDXXX  <- HMDXXX[!is.na(HMDXXX[,"ST_Per_LE_FY"]), ]
	HMDXXX$Subpop.Code
}

############################################################################
# 3) getJMDprefectures()
############################################################################

#' @title get a named vector of JMD prefecture codes
#' 
#' @description This is a helper function for those familiar with prefecture names but not with prefecture codes (and vice versa). It is also useful for looped downloading of data.
#' 
#' @return a character vector of 2-digit prefecture codes. Names correspond to the proper names given in the English version of the HMD webpage.
#' 
#' @importFrom XML readHTMLTable
#' 
#' @export
#' 
#' @examples \dontrun{ (prefectures <- getJMDprefectures()) }
#' 
getJMDprefectures <- function(){
	Prefs <- as.matrix(XML::readHTMLTable("http://www.ipss.go.jp/p-toukei/JMD/index-en.html",
					which = 1, stringsAsFactors = FALSE, skip.rows = c(1:4)))
	# get codes. rows read from left to right
	Prefectures <- sprintf("%.2d", row(Prefs) * 4 + col(Prefs) - 5)
	names(Prefectures) <- c(Prefs)
	Prefectures[order(Prefectures)]
}

############################################################################
# 4) getCHMDprovinces()
############################################################################

#' @title get a named vector of CHMD province codes
#' 
#' @description This is a helper function to get a vector of 3-character province codes.
#' 
#' @return a character vector of 3 character province codes.
#' 
#' @export
#' 
#' @examples \dontrun{ (provs <- getCHMDprovinces()) }
#' 
getCHMDprovinces <- function(){
	# it's a small list, so why both scraping?-- include "can" for posterity.
	sort(c("can","nfl","pei","nsc","nbr","que","ont","man","sas","alb","bco","nwt","yuk"))
}
