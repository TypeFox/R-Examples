#' Validate Finnish personal identification numbers (hetu) 
#'
#' @param hetu Finnish personal identification number as a character vector, 
#' 	  or vector of identification numbers as a character vectors.
#'
#' @return Is the given string a valid Finnish personal identification number, 
#' 	   \code{TRUE} or \code{FALSE}.
#' 
#' @author Jussi Paananen \email{louhos@@googlegroups.com}
#' 
#' @seealso \code{\link{hetu}} For extracting information from Finnish 
#' 	    personal identification numbers. 
#' 
#' @examples
#' valid_hetu("010101-0101") # TRUE
#' valid_hetu("010101-010A") # FALSE
#' @export

valid_hetu <- function(hetu) {
  # Try to create hetu-object from the given hetu, check if created object 
  # is of the correct class 
  if (length(hetu) > 1) {
    return(sapply(hetu, FUN=valid_hetu))
  }

  return(class(hetu(hetu)) == "data.frame")
}


#' Extract information from Finnish personal identification numbers (hetu)
#'
#' @param hetu Finnish personal identification number as a character vector, 
#' 	  or vector of identification numbers as a character vectors
#' @param extract Extract only selected part of the information. 
#'    Valid values are "\code{hetu}", "\code{gender}", "\code{personal.number}",
#'    "\code{checksum}", "\code{date}", "\code{day}", "\code{month}", 
#'    "\code{year}", "\code{century.char}".
#'    If \code{NULL} (default), returns all information. 
#'
#' @return Finnish personal identification number data.frame,
#'         or if extract parameter is set, the requested part of the 
#'	   information as a vector. Returns \code{NA} if the given character 
#'	   vector is not a valid Finnish personal identification number.
#' \item{hetu}{Finnish personal identification number as a character vector.}
#' \item{gender}{Gender of the person as a character vector ("Male" or "Female").}
#' \item{personal.number}{Personal number part of the identification number.}
#' \item{checksum}{Checksum for the personal identification number.}
#' \item{date}{Birthdate.}
#' \item{day}{Day of the birthdate.}
#' \item{month}{Month of the birthdate.}
#' \item{year}{Year of the birthdate.}
#' \item{century.char}{Century character of the birthdate: 
#'                     + (1800), - (1900) or A (2000). }
#' 
#' @author Jussi Paananen \email{louhos@@googlegroups.com}
#' 
#' @seealso \code{\link{valid_hetu}} For validating Finnish personal 
#' 	    identification numbers.
#' @examples
#' hetu("111111-111C")
#' hetu("111111-111C")$date
#' hetu("111111-111C")$gender
#' # Same as previous, but using extract argument
#' hetu("111111-111C", extract="gender")
#' 
#' # Process a vector of hetu's
#' hetu(c("010101-0101", "111111-111C"))
#' 
#' # Process a vector of hetu's and extract gender information from each
#' hetu(c("010101-0101", "111111-111C"), extract="gender")
#' @export

hetu <- function(hetu, extract=NULL) {
  
  if (!is.null(extract)) {
    if (!extract %in% c("hetu", "gender", "personal.number", "checksum", 
       		        "date", "day", "month", "year", "century.char")) {
      stop("Trying to extract invalid part of hetu")
    }
  }
  
  # Check if the input parameter is a vector
  if (length(hetu) > 1) {
    if (is.null(extract)) {
      res <- lapply(hetu, FUN=hetu, extract=extract)
      # Convert dates to characters to avoid conversion problems
      for (i in 1:length(res)) {res[[i]]$date <- as.character(res[[i]]$date)}
      # Convert list to data.frame
      res <- do.call(rbind.data.frame, res) 
      # dates back to dates
      res$date <- as.Date(as.character(res$date))
      # Return
      return(res)
    } else {
      return(unname(do.call("c", lapply(hetu, FUN=hetu, extract=extract))))
    }    
  }
  
  # Convert to character vector if necessary
  if(!is.character(hetu)) hetu <- as.character(hetu)
  
  # Check general format
  match <- regexpr("^[0-9]{6}[\\+-A][0-9]{3}[0123456789ABCDEFHJKLMNPRSTUVWXY]$", hetu)
  if (match == -1 ) {
    return(NA)
  }
  
  # Check day
  day <- as.numeric(substr(hetu, start=1, stop=2))
  if (!((day >= 1) && (day <= 31))) {
    return(NA)
  }
  
  # Check month
  month <- as.numeric(substr(hetu, start=3, stop=4))
  if (!((month >= 1) && (month <= 12))) {
    return(NA)
  }
  
  # Check year
  year <- as.numeric(substr(hetu, start=5, stop=6))
  if (!((year >= 00) && (year <= 99))) {
    return(NA)
  }
  
  # Check century
  century <- substr(hetu, start=7, stop=7)
  if (!century %in% c("+", "-", "A")) {
    return(NA)
  }
  
  # Construct complete year from century character and 2-digit year
  
  ## Pad leading zero to a 2-digit year if needed
  year <- formatC(year, flag=0, width=2) 
  
  if (century == "+") {
    full.year <- as.numeric(paste("18", year, sep=""))
  }
  if (century == "-") {
    full.year <- as.numeric(paste("19", year, sep=""))
  }
  if (century == "A") {
    full.year <- as.numeric(paste("20", year, sep=""))
  }
  
  # Check if date exists
  date <- as.Date(paste(day, "/", month, "/", full.year, sep=""), "%d/%m/%Y")
  if (is.na(date)) {
    return(NA)
  }
  
  # Check personal identification number
  personal <- as.numeric(substr(hetu, start=8, stop=10))
  if (!((personal >= 2) && (personal <= 899))) {
    return(NA)
  }
  
  # Check checksum character validity
  check <- substr(hetu, start=11, stop=11)
  checklist <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F", "H", "J", "K", "L", "M", "N", "P", "R", "S", "T", "U", "V", "W", "X", "Y")
  names(checklist) <- 0:30
  if (!check %in% checklist) {
    return(NA)
  }
  # Check checksum character
  mod <- as.numeric(paste(substr(hetu, start=1, stop=6), 
      	 		substr(hetu, start=8, stop=10), sep="")) %% 31
  if (check != checklist[as.character(mod)]) {
    return(NA)
  }
  
  # Check gender
  if ((personal %% 2) == 0) {
    gender <- "Female"
  } else {
    gender <- "Male"
  }

  # Create hetu-object
  object <- list(hetu = hetu, gender=gender, 
  	         personal.number=personal, 
  	         checksum=check, date=date, day=day, month=month, 
		 year=full.year, century.char=century)
  #class(object) <- "hetu"
  
  # Return full object or only requested part
  if (is.null(extract)) {
    return (as.data.frame(object))
  }
  else {
    return(unname(do.call("c", object[extract])))
  }
}

