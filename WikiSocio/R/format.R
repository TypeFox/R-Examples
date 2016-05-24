divide_list <- function(x, number) {
    if (length(x) > number) {
        rest <- length(x)%%number
        if (rest == 0) {
            result <- matrix(x, ncol = number, byrow = TRUE)
            result <- split(result, row(result))
        } else {
            index <- length(x) - rest
            result <- matrix(x[1:index], ncol = number, byrow = TRUE)
            result <- split(result, row(result))
            result[[length(result) + 1]] <- x[index + 1:length(x)]
            length(result[[length(result)]]) <- length(x) - index
        }
        return(result)
    } else {
        return(as.vector(x))
    }
}

#' Query a page by xpath means
#'
#' @param url URL of the page. Please note that this is different of the title of the wiki page.
#' @param xpath XPATH query to execute
#'
#' @import XML
#' @import RCurl
#' 
#' @return A character vector
#' @export

xpath <- function(url, xpath) {
    
    dom <- htmlParse(getURL(url, followlocation = TRUE, httpheader = c(`User-Agent` = "WikiSocio"), encoding = "UTF-8"))
    return(xpathSApply(dom, xpath, xmlValue))
    
}

#' Recode a date vector into periods
#'
#' @param date A timestamp vector
#' @param period The breaks witch \code{date} sould be cut with
#'
#' @return a numeric vector with the periods replacing the date
#' @export
#'
period <- function(date, period = c("2005-05-01", "2009-12-01")) {
    
    if (period == "year") {
        period <- c("2001-01-01", "2002-01-01", "2003-01-01", "2004-01-01", "2005-01-01", "2006-01-01", "2007-01-01", "2008-01-01", "2009-01-01", 
            "2010-01-01", "2011-01-01", "2012-01-01", "2013-01-01", "2014-01-01", "2015-01-01", "2016-01-01")
    }
    
    period <- tryCatch(as.Date(period), error = function(e) stop("period is not a time vector"))
    date <- tryCatch(as.Date(date), error = function(e) NA)
    
    return(tryCatch(cut(date, breaks = c(as.Date("1900-01-01"), period, as.Date("3000-01-01"))), error = function(e) NULL))
} 
