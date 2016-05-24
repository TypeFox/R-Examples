readHMDformat <-
function(CNTRY = NULL, username = NULL, password = NULL, fixup = TRUE){
    age2int <- function(Age){
        ## replaces + and - with nothing
        ## e.g., "99+" becomes "99"
        as.integer(gsub("[-+]", "", unlist(lapply(strsplit(as.character(Age),split="-"),"[[",1))))
    }


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
    userInput <- function(silent = FALSE){
        if (interactive()){
            if(!silent){
                cat("\ntype in a single character string.")
            }
            out <- scan(file = "", n = 1, what = "character")
        } else {
            stop("User input only works for an interactive R session")
        }
        invisible(out)
    }

    

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
    dths <- readHMDweb(CNTRY=CNTRY, item="Deaths_1x1", username=username, password=password)[,1:5]
    expos <- readHMDweb(CNTRY=CNTRY, item="Exposures_1x1", username=username, password=password)[,1:5]
    return(list(deaths=dths, exposures=expos))
}
