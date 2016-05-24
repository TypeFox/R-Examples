.onLoad <- function(lib, pkg) {
	VariABEL.version <- "0.9-1"
#	packageStartupMessage("VariABEL v.",VariABEL.version,"(November 07, 2013) loaded\n")
	
	# check for updates and news
	address <- c(
			"http://genabel.r-forge.r-project.org/version_and_news.html",
			"http://www.genabel.org/sites/default/files/version_and_news.html"
	)
	svtmo <- options("timeout")
	options("timeout"=10)
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
			a <- tolower(fulltext)
			a <- a[grep("<vastable>",a)+1]
			if (length(a)>0) {
				# message to all users
				strnews <- grep("<messagetoall>",tolower(fulltext))
				endnews <- grep("</messagetoall>",tolower(fulltext))
				if (length(strnews)>0 && length(endnews)>0) 
					if ((endnews-1) >= (strnews+1)) {
#						packageStartupMessage(fulltext[(strnews+1):(endnews-1)],sep="\n")
					}
				# compare versions
				a <- strsplit(a,"")[[1]]
				ver <- a[grep("[0-9]",a)]
				ver <- paste(ver[1],".",ver[2],"-",ver[3],sep="")
				if (VariABEL.version != ver) {
#					packageStartupMessage(  "\nInstalled VariABEL version (",VariABEL.version,") is not the same as stable\n",
#							"version available from CRAN (",ver,"). Unless used intentionally,\n",
#							"consider updating to the latest CRAN version. For that, use\n",
#							"'install.packages(\"VariABEL\")', or ask your system administrator\n",
#							"to update the package.\n\n",sep="")
					# check for new-version news
					strnews <- grep("<vanews>",tolower(fulltext))
					endnews <- grep("</vanews>",tolower(fulltext))
					if (length(strnews)>0 && length(endnews)>0) 
						if ((endnews-1) >= (strnews+1)) {
#							packageStartupMessage(fulltext[(strnews+1):(endnews-1)],sep="\n")
						}
				}
			}
		}
		#rm(a,fulltext,ver)
	}
	options("timeout"=svtmo)
	#rm(tryRes0,tryRes1,conn,svtmo)
}
