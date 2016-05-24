#' Getting the list of the contribution of a contributor.
#'
#' @param x The name of the contributor
#' @param namespace The namespace of contributions
#' @param domain The domain of the wiki
#' @param page If \code{TRUE}, contrib_list return only a list of all different page in witch the contribors made a modification. If \code{FALSE}, return a data-frame with all the contributions, their timestamp and their weight
#'
#' @importFrom httr GET content
#'
#' @return Depending the value of page, weither a character vector containing the names of all contributors, or a data-frame containing all the revisions with the name of the contributor, a timestamp and a weight
#' @export
#' 
#'
#' @examples
#' 
#' # All the contribution of an user of the french wiki.
#' contrib_list('cafeine05') 
#' 
#' # Return a character vector with all the page modified by this contributor.
#' contrib_list('cafeine05',page=TRUE) 
contrib_list <- function(x, namespace = "0", domain = "fr", page = F) {
    uccontinue <- NULL
    result <- data.frame()
    unname(result)
    
    repeat {
        if (is.null(uccontinue)) {
            query = list(action = "query", list = "usercontribs", format = "json", ucnamespace = namespace, ucprop = "title|timestamp|sizediff", 
                uclimit = "max", ucuser = x)
        } else {
            query = list(action = "query", list = "usercontribs", format = "json", ucnamespace = namespace, ucprop = "title|timestamp|sizediff", 
                uclimit = "max", ucuser = x, uccontinue = uccontinue)
        }
        
        exec <- GET(paste("https://", domain, ".wikipedia.org/w/api.php", sep = ""), query = query)
        
        content <- content(exec, "parsed")
        uccontinue <- tryCatch(content$continue$uccontinue, error = function(e) NULL)
        dF <- tryCatch(content$query$usercontribs, error = function(e) NULL)
        
        if (!is.null(dF)) {
            
            dF <- lapply(dF, function(x) {
                if (!"sizediff" %in% names(x)) {
                  x$sizediff <- NA
                }
                if (!"timestamp" %in% names(x)) {
                  x$timestamp <- NA
                }
                if (!"title" %in% names(x)) {
                  x$title <- NA
                }
                matrix(c(x$title, x$timestamp, x$sizediff), ncol = 3, byrow = FALSE)
            })
            
            dF <- do.call(rbind, dF)
            
            dF <- as.data.frame(dF)
            unname(dF)
            
            result <- rbind(result, dF)
            
        }
        
        if (is.null(uccontinue)) {
            break
        }
    }
    
    if (page) {
        result <- as.vector(unique(result[, 1]))
    }
    return(result)
    
} 
