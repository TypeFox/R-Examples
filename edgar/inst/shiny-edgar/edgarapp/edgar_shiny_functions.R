# function to download file and return FALSE if download error
DownloadFile <- function(link, dfile) {
	tryCatch({
		utils::download.file(link, dfile, quiet = TRUE)
		return(TRUE)
		}, error = function(e) {
				return(FALSE)
	})
}
	
# Master Index Download function
GetMasterIndexShiny <- function(year.array) {
    status.array <- data.frame()
    
    for (i in 1:length(year.array)) {
        year <- year.array[i]
        year.master <- data.frame()
        quarterloop <- 4
        if (year == format(Sys.Date(), "%Y")) {
            quarterloop <- ceiling(as.integer(format(Sys.Date(), "%m"))/3)
        }
        
        for (quarter in 1:quarterloop) {
            # save downloaded file with specific name
            dfile <- paste0("Master Index/", year, "QTR", quarter, "master.gz")
            file <- paste0("Master Index/", year, "QTR", quarter, "master")
            # form a link to download master file
            link <- paste0("ftp://ftp.sec.gov/edgar/full-index/", year, "/QTR", quarter, "/master.gz")
            # progress bar
            incProgress(1/quarter, detail = paste0("Yr-", year, ", Qtr-", quarter))
            
            res <- DownloadFile(link, dfile)
            if (res) {
                # Unzip gz file
                R.utils::gunzip(dfile, destname = file, temporary = FALSE, skip = FALSE, overwrite = TRUE, remove = TRUE)
                cat("Successfully downloaded Quarter Master File", file, " for year:", year, "and quarter:", quarter, "...")
                
                # Removing "'" so that scan with "|" not fail due to occurrence of "'" in company name
                data <- gsub("'", "", readLines(file))
                # writting back to storage
                writeLines(data, file)
                
                d <- scan(file, what = list("", "", "", "", ""), flush = F, skip = 10, sep = "|")
                data <- data.frame(CIK = d[[1]], COMPANY_NAME = d[[2]], FORM_TYPE = d[[3]], DATE_FILED = d[[4]], EDGAR_LINK = d[[5]], QUARTER = quarter)
                year.master <- rbind(year.master, data)
                file.remove(file)
                status.array <- rbind(status.array, data.frame(Filename = paste0(year, ": quarter-", quarter), status = "Download success"))
            } else {
                status.array <- rbind(status.array, data.frame(Filename = paste0(year, ": quarter-", quarter), status = "Server Error"))
            }
        }
        assign(paste0(year, "master"), year.master)
        save(year.master, file = paste0("Master Index/", year, "master.Rda"))
    }
    return(status.array)
}


# Function for forming word frequency
GetwordfrqShiny <- function(filpath) {
    text <- readLines(filpath)
    text <- paste(text, collapse = " ")
    
    # Extract text from html file
    doc <- XML::htmlParse(text, asText = TRUE)
    text <- XML::xpathSApply(doc, "//text()[not(ancestor::script)][not(ancestor::style)][not(ancestor::noscript)][not(ancestor::form)]", 
							 XML::xmlValue)
    text <- paste(text, collapse = " ")
    
    # convert into corpus
    text <- tm::Corpus(tm::VectorSource(text))
    # clean text
    cleantext <- function(data.text.corpus) {
        data.text.corpus <- tm::tm_map(data.text.corpus, tm::removePunctuation)  # Remove punctuation marks
        data.text.corpus <- tm::tm_map(data.text.corpus, tm::removeNumbers)  # Remove punctuation marks
        data.text.corpus <- tm::tm_map(data.text.corpus, tm::stripWhitespace)  # Remove Numbers
        data.text.corpus <- tm::tm_map(data.text.corpus, function(x) tm::removeWords(x, tm::stopwords()))  # Remove stop words
        data.text.corpus <- tm::tm_map(data.text.corpus, tm::content_transformer(tolower))  # Convert text to lower case
        return(data.text.corpus)
    }
    text <- cleantext(text)
    word.frq <- tm::termFreq(text[[1]])
    word.frq <- data.frame(WORD = names(word.frq), FREQUENCY = word.frq, row.names = NULL)
    word.frq <- word.frq[order(-word.frq$FREQUENCY), ]
    rownames(word.frq) <- NULL
    return(word.frq)
}

# Function to create links 
CreateLinkShiny <- function(d) {
    links <- (d$EDGAR_LINK)
    links <- gsub(".txt", "-index.htm", links)
    d$SELECT <- sprintf("<a href=\"%s\">%s</a>", paste0("http://www.sec.gov/Archives/", links), "VIEW")
    d$EDGAR_LINK <- NULL  # Remove edgar links 
    return(d)
}

# Function for downloading Filings
DownloadFilingsShiny <- function(year, cik.no, form.type) {
    yr.master <- paste0(year, "master.Rda")
    stat.Filing2 <- data.frame()
    
    # function to download file and return FALSE if download error
    DownlFilingshiny <- function(LINK, dest.filename) {
        tryCatch({
            utils::download.file(LINK, dest.filename, quiet = TRUE)
            return(TRUE)
        }, error = function(e) {
            return(FALSE)
        })
    }
    
    if (file.exists(paste0("Master Index/", yr.master))) {
        load(paste0("Master Index/", yr.master))
        
        if (nrow(year.master) > 0) {
            # check if user want to download all cik or specific cik
            if (!cik.no == "ALL") {
                year.master <- year.master[year.master$CIK == cik.no, ]
            }
            # if cik.no not found in the master file then show msg and exit
            if (nrow(year.master) == 0) {
                msg3 <- paste0("CIK No: ", cik.no, " not found in the selected year")
                err <- tcltk::tkmessageBox(message = msg3, icon = "error")
                # stop(msg3)
                return(stat.Filing2)
            }
            
            if (!form.type == "ALL") {
                year.master <- year.master[year.master$FORM_TYPE == form.type, ]
            }
            
            # if form_type not found in the master file then show msg and exit
            if (nrow(year.master) == 0) {
                msg2 <- paste0("Form Type: ", form.type, " not Filed by selected CIK")
                err <- tcltk::tkmessageBox(message = msg2, icon = "error")
                # stop(msg2)
                return(stat.Filing2)
            }
            
            # downloading files
            total.files <- nrow(year.master)
            msg3 <- paste0("Total EDGARs to be downloaded:", total.files, "\nDo you want to download the EDGAR files")
            choice <- tcltk::tkmessageBox(message = msg3, type = "yesno", default = "no")
            
            if (as.character(choice) == "yes") {
                # progress bar
                dir.create("Edgar filings")
                progress.bar <- tcltk::tkProgressBar(title = "Progress Bar", min = 0, max = total.files, width = 400)
                f.type <- gsub("/", "", form.type)
                new.dir <- paste0("Edgar filings/", cik.no, "_", f.type, "_", year)
                dir.create(new.dir)
                
                for (i in 1:total.files) {
                  LINK <- paste0("http://edgar.sec.gov/Archives/", year.master$EDGAR_LINK[i])
                  f.type <- gsub("/", "", year.master$FORM_TYPE[i])
                  dest.filename <- paste0(new.dir, "/", year.master$CIK[i], "_", 
				                          f.type, "_", year.master$DATE_FILED[i], ".txt")
                  cat("Downloading file:", LINK, "\n")
                  res <- DownlFilingshiny(LINK, dest.filename)
                  
                  if (res) {
                    temp.status <- data.frame(Link = LINK, Status = "Download success")
                  } else {
                    temp.status <- data.frame(Link = LINK, Status = "Server Error")
                  }
                  
                  stat.Filing2 <- rbind(stat.Filing2, temp.status)
                  prg <- paste0(ceiling(i/total.files * 100), "% done")
                  
                  tcltk::setTkProgressBar(progress.bar, i, label = prg)
                }
                close(progress.bar)
                return(stat.Filing2)
                
            }
        } else {
            msg3 <- "Rda file is corrupted. Please redownload the master Index file for the selected year using 'Get Master Index' tab."
            err <- tcltk::tkmessageBox(message = msg3, icon = "error")
            return(stat.Filing2)
        }
    } else {
        errmsg <- paste0("Current directory does not contains ", yr.master, " file in 'Master Index' directory. Please download Master files using 'Get Master Index tab'.")
        err <- tcltk::tkmessageBox(message = errmsg, icon = "error")
        return(stat.Filing2)
    }
}

# function for downloading daily Index
GetDailyInfoShiny <- function(day, month, year) {
    date <- paste0(year, month, day)
    filename <- paste0("Daily Index/daily_idx_", date)
    link1 <- paste0("ftp://ftp.sec.gov/edgar/daily-index/master.", date, ".idx")
    link2 <- paste0("ftp://ftp.sec.gov/edgar/daily-index/", year, "/QTR", 
	                ceiling(as.integer(month)/3), "/master.", date, ".idx.gz")
    link3 <- paste0("ftp://ftp.sec.gov/edgar/daily-index/", year, "/QTR", 
	                ceiling(as.integer(month)/3), "/master.", substr(as.character(year), 3, 4), month, day, ".idx")
    link4 <- paste0("ftp://ftp.sec.gov/edgar/daily-index/", year, "/QTR", 
	                 ceiling(as.integer(month)/3), "/master.", date, ".idx")
    down.success = FALSE
    
    if (year < 1999) {
        fun.return3 <- DownloadFile(link3, filename)
        if (fun.return3 == 1) {
            down.success = TRUE
        }
    }
    
    if (year > 1998 && year < 2012) {
        fun.return4 <- DownloadFile(link4, filename)
        if (fun.return4 == 1) {
            down.success = TRUE
        }
    }
    
    if (year > 2011) {
        fun.return1 <- DownloadFile(link1, filename)
        if (fun.return1 == 1) {
            down.success = TRUE
        } else {
            fun.return2 <- DownloadFile(link2, filename)
            if (fun.return2 == 1) {
                down.success = TRUE
            }
        }
    }
    
    if (down.success) {
	
	    # Removing "'" so that scan with "|" not fail due to occurrence of "'" in company name
        temp.data <- gsub("'", "", readLines(filename))
        # writting back to storage
        writeLines(temp.data, filename)
		
        d <- scan(filename, what = list("", "", "", "", ""), flush = F, skip = 10, sep = "|")
        data <- data.frame(CIK = d[[1]], COMPANY_NAME = d[[2]], FORM_TYPE = d[[3]], 
						   DATE_FILED = paste0(year, "-", month, "-", day), EDGAR_LINK = d[[5]])
        data <- CreateLinkShiny(data)
        return(data)
    } else {
        return(0)
    }
}
