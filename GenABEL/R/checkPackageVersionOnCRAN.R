#' checks what is the version of package on CRAN
#' 
#' Checks what is the version of package on CRAN.
#' The CRAN page (baseUrlCRAN+packageName) is checked 
#' and parsed extracting the line with
#' "Package source:	 packageName_Version.tar.gz" 
#' e.g. 
#' "Package source:	 GenABEL_1.6-9.tar.gz"
#' and then the 'Version' is returned. 
#' Otherwise, NULL is returned. 
#' 
#' @return string containing CRAN version 
#' of the package
#' 
#' @param packageName name of the package to check
#' @param baseUrlCRAN path to CRAN repository
#' @param timeout web chack timeout
#' 
#' @examples 
#' library(GenABEL)
#' packageVersion("GenABEL")
#' checkPackageVersionOnCRAN("GenABEL")
#' 
#' @author Yurii Aulchenko
#'
checkPackageVersionOnCRAN <- function(packageName,baseUrlCRAN="http://cran.r-project.org/web/packages/", 
		timeout = 2)
{
	# change default timout
	svtmo <- options("timeout")
	options("timeout"=timeout)
	# page to check is
	pageAddress <- paste(baseUrlCRAN,packageName,sep="/")
	# establish connection to the CRAN page of the package
	suppressWarnings(
			conn <- try( url(pageAddress) , silent=TRUE )
	)
	# if connection ok, read full page, store the results in pageContent; if failed, pageContent <- "try-error"
	if ( all( class(conn) != "try-error") ) {
		suppressWarnings(
				pageContent <- try( readLines(conn) , silent=TRUE )
		)
		close(conn)
	} else {
		pageContent <- "try-error"
		class(pageContent) <- "try-error"
	}
	# restore default timeout
	options("timeout"=svtmo)
	# if failed in reading (pageContent is "try-error"), return NULL
	if (class(pageContent) == "try-error") return(NULL)
	# parse the page and get string starting with "Package source:"
	targetLine <- pageContent[grep("source:",pageContent)]
	# split the string at "Package_" and ".tar.gz"; the element before the last will contain the version
	splitPattern <- paste(packageName,"_|.tar.gz",sep="")
	stringSplit <- strsplit(targetLine,splitPattern)
	cranVersion <- stringSplit[[1]][length(stringSplit[[1]])-1]
	# return version
	return(cranVersion)
}
