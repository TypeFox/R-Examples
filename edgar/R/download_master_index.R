#' Retrieves quarterly Master Index.
#'
#' \code{DownloadMasterIndex} retrieves the quarterly Master Index from SEC site.
#'
#' DownloadMasterIndex function takes 'year' as an input parameter from user,  
#' asks the user to locate working directory, download quarterly master index
#' from http://edgar.sec.gov/Archives, strips the headers, converts it into tabular 
#' form, and merges such quarterly tables into yearly tables in RData format.
#' Function creates new directory 'Master Index' into working directory 
#' to save these Rda Master Index. Please note, for all other functions in this 
#' package needs to locate the same working directory to access these Rda index files.  
#'  
#' @param year.array year in integer or integer array containing years for which Master 
#' Index are to be downloaded.
#' 
#' @return Function retrieves quarterly Master Index files 
#' from \url{http://edgar.sec.gov/Archives} site and returns download status dataframe.
#'   
#' @examples
#' \dontrun{
#' 
#' report <- DownloadMasterIndex(1995) 
#' ## Download quarterly Master Index files for the year 1990 and stores into yearly  
#' ## 1995master.Rda file. It returns download report in dataframe format.
#' 
#' report <- DownloadMasterIndex(c(1994, 1995, 2006)) 
#' ## Download quarterly Master Index files for the years 1994, 1995, 2006 and stores into 
#' ## different {year}master.Rda files. It returns download report in dataframe format.
#'}

DownloadMasterIndex <- function(year.array) {
    if (!is.numeric(year.array)) {
        msg <- "Please enter valid year"
        err <- tcltk::tkmessageBox(message = msg, icon = "error")
        stop(msg)
    }
    
    # function to download file and return FALSE if download error
    DownloadFile <- function(link, dfile) {
        tryCatch({
            utils::download.file(link, dfile, quiet = TRUE)
            return(TRUE)
        }, error = function(e) {
            return(FALSE)
        })
    }
    
    options(warn = -1)
    cat("Select working directory in pop up menu\n")
    setwd(jchoose.dir(default = getwd(), caption = "Choose working directory"))
    dir.create("Master Index")
    
    status.array <- data.frame()
    for (i in 1:length(year.array)) {
        year <- year.array[i]
        year.master <- data.frame()
        quarterloop <- 4
        if (year == format(Sys.Date(), "%Y")) {
            quarterloop <- ceiling(as.integer(format(Sys.Date(), "%m"))/3)
        }
        
        for (quarter in 1:quarterloop) {
            # save downloaded file as specific name
            dfile <- paste0("Master Index/", year, "QTR", quarter, "master.gz")
            file <- paste0("Master Index/", year, "QTR", quarter, "master")
            
            # form a link to download master file
            link <- paste0("ftp://ftp.sec.gov/edgar/full-index/", year, "/QTR", quarter, "/master.gz")
            
            res <- DownloadFile(link, dfile)
            if (res) {
                # Unzip gz file
                R.utils::gunzip(dfile, destname = file, temporary = FALSE, skip = FALSE, overwrite = TRUE, remove = TRUE)
                cat("Successfully downloaded Master Index for year:", year, "and quarter:", quarter)
                
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
