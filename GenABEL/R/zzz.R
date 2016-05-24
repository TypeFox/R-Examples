DISABLED.onAttach <- function(lib, pkg) {
	## this is something which should be fixed: both version 
	## and date can come from DESCRIPTION!
	#pkgDescription <- packageDescription(pkg)
	#pkgVersion <- pkgDescription$Version
	#pkgDate <- pkgDescription$Date
	pkgVersion <- "1.7-7"
	pkgDate <- "June 29, 2013"
	welcomeMessage <- paste(pkg," v. ",pkgVersion," (",pkgDate,") loaded\n",sep="")
	# check if CRAN version is the same as loaded
	cranVersion <- try( checkPackageVersionOnCRAN(pkg) )
	if (!is.null(cranVersion) & !( class(cranVersion) == "try-error") ) 
		if (pkgVersion != cranVersion) {
			welcomeMessage <- paste(welcomeMessage,
					"\nInstalled ",pkg," version (",pkgVersion,") is not the same as stable\n",
					"version available from CRAN (",cranVersion,"). Unless used intentionally,\n",
					"consider updating to the latest CRAN version. For that, use\n",
					"'install.packages(\"",pkg,"\")', or ask your system administrator\n",
					"to update the package.\n\n",sep="")
		}
	# check for news
	address <- c(
			"http://genabel.r-forge.r-project.org/version_and_news.html",
			"http://www.genabel.org/sites/default/files/version_and_news.html"
	)
	svtmo <- options("timeout")
	options("timeout"=2)
	tryRes1 <- 0; class(tryRes1) <- "try-error"
	curaddr <- 1
	while (class(tryRes1) == "try-error" && curaddr <= length(address) ) {
		suppressWarnings(
				tryRes0 <- try(conn <- url(address[curaddr]),silent=TRUE)
		)
		suppressWarnings(
				tryRes1 <- try(fulltext <- readLines(conn),silent=TRUE)
		)
		close(conn)
		curaddr <- curaddr + 1
	}
	if (class(tryRes1) != "try-error") {
		if (length(fulltext)>0)
		{
			# message to all users
			strnews <- grep("<messagetoall>",tolower(fulltext))
			endnews <- grep("</messagetoall>",tolower(fulltext))
			if (length(strnews)>0 && length(endnews)>0) 
				if ((endnews-1) >= (strnews+1)) {
					welcomeMessage <- paste(welcomeMessage,
							fulltext[(strnews+1):(endnews-1)],sep="\n")
				}
			# check for specific package news
			strnews <- grep(paste("<",pkg,"news>",sep=""),tolower(fulltext))
			endnews <- grep(paste("</",pkg,"news>",sep=""),tolower(fulltext))
			if (length(strnews)>0 && length(endnews)>0) 
				if ((endnews-1) >= (strnews+1)) {
					welcomeMessage <- paste(welcomeMessage,
							fulltext[(strnews+1):(endnews-1)],sep="\n")
				}
		}
		#rm(a,fulltext,ver)
	}
	options("timeout"=svtmo)
	packageStartupMessage(welcomeMessage)
}
