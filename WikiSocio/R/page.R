#' List of inter-articles links on a wiki page
#'
#' @param x The title of the page
#' @param domain The domain of the wiki where the page is located
#' @param namespace The namespace where the page need to be to be ept in the graph
#'
#' @return A character vector
#' @export
#' 
#' @import pbapply
#' 
#' @family page functions
#'
#' @examples
#' page_links('Action') # listing all the links contained in the 'action' article.
page_links <- function(x, domain = "fr", namespace = "0") {
    
    plcontinue <- NULL
    result <- vector(mode = "character")
    repeat {
        
        if (is.null(plcontinue)) {
            query = list(action = "query", prop = "links", format = "json", plnamespace = namespace, pllimit = "max", titles = x, redirects = "")
        } else {
            query = list(action = "query", prop = "links", format = "json", plnamespace = namespace, pllimit = "max", titles = x, plcontinue = plcontinue, 
                redirects = "")
        }
        
        exec <- GET(paste("https://", domain, ".wikipedia.org/w/api.php", sep = ""), query = query)
        
        content <- content(exec, "parsed")
        plcontinue <- content$continue$plcontinue
        data <- content$query$pages
        idPage <- data[[1]][1]
        data <- data[[as.character(idPage)]][["links"]]
        data <- sapply(data, function(x) {
            x[[2]]
        })
        result <- c(result, data)
        
        if (is.null(plcontinue)) {
            break
        }
    }
    
    return(list(list = result, nb = length(result)))
}

#' Test if page(s) is(are) linked
#'
#' @param from The title of the page from the link is supposed to go from
#' @param to Either a string or a character vector representing the page(s) where the link(s) is(are) supposed to go to.
#' @param namespace The namespace of the page to look.
#' @param domain The domain where the wiki is located
#'
#' @return A character vector.
#' @export
#' 
#' @family page functions
#'
#' @examples
#' # Testing if the article 'Action' contains a link to 'Sociologie' in the french wiki.
#' page_islink('Action','Sociologie') 

page_islink <- function(from, to, namespace = "0", domain = "fr") {
    
    
    to <- divide_list(to, 50)
    to <- lapply(to, paste, collapse = "|")
    
    result <- lapply(to, function(x) {
        data <- GET(paste("https://", domain, ".wikipedia.org/w/api.php", sep = ""), query = list(action = "query", prop = "links", format = "json", 
            pllimit = "max", pltitles = x, titles = from, redirects = ""))
        
        data <- content(data, "parsed")$query$pages
        idPage <- data[[1]][1]
        data <- data[[as.character(idPage)]]["links"]
        
        if (!is.null(unlist(data))) {
            data <- matrix(unlist(data), ncol = 2, byrow = TRUE)
            data[, 2]
        }
    })
    return(as.vector(unlist(result)))
}

#' List all categories where the page is located
#'
#' @param title The title of the page
#' @param domain The domain where the wiki is located
#'
#' @return A character vector giving 
#' @export
#'
#' @examples
#' # Downloading the list of categories in the action page
#' page_cat('Action') 
page_cat <- function(title, domain = "fr") {
    clcontinue <- NULL
    result <- vector(mode = "character")
    repeat {
        
        if (is.null(clcontinue)) {
            query = list(action = "query", prop = "categories", format = "json", cllimit = "max", titles = title, redirects = "")
        } else {
            query = list(action = "query", prop = "categories", format = "json", cllimit = "max", titles = title, clcontinue = clcontinue, redirects = "")
        }
        
        exec <- GET(paste("https://", domain, ".wikipedia.org/w/api.php", sep = ""), query = query)
        parsed <- content(exec, "parsed")
        
        content <- parsed$query$pages[[1]]$categories
        
        if (length(content) > 0) {
            
            content <- matrix(unlist(content), ncol = 2, byrow = TRUE)[, 2]
            clcontinue <- parsed$continue$clcontinue
            result <- c(result, content)
            
        } else {
            
            return(NULL)
        }
        
        if (is.null(clcontinue)) {
            break
        }
    }
    return(result)
}

#' Giving the proportion of anonmous contribution in a revisions list
#'
#' @param revisions Revision list built with page_revisions
#'
#' @return An integer.
#' @export
#' 
#' @family page functions
#'
#' @examples
#' # Downloading the list of the revisions of the article 'Action' in the french wiki
#' revisions <- page_revisions('Action') 
#' 
#' # Returning the percentage of anonymous contribution in this article
#' page_anon(revisions) 
page_anon <- function(revisions) {
    revisions <- as.data.frame(revisions)
    return((length(revisions[revisions[, 4] == TRUE, 4])/length(revisions[, 4])) * 100)
}

#' Downloading the list of contributions for one page
#'
#' @param page The title of the page
#' @param domain The domain where the wiki is located
#' @param contrib A boolean indicating whether to return or not only contributors' names
#'
#' @return A data-frame containing the username of the user (or the IP if anonymous contribution), the timestamp, the size of the revision, a boolean indicating weither the contribution is anonymous or not, and the difference beetween the contribution and the previous
#' @export
#'
#' @family page functions
#'
#' @examples
#' # Downloading the list of contribution for the 'action' page in the french wiki
#' page_revisions('Action') 
page_revisions <- function(page, domain = "fr", contrib = F) {
    continue <- NULL
    result <- data.frame(matrix(ncol = 5))
    names(result) <- c("user", "timestamp", "size", "anon", "weight")
    repeat {
        if (is.null(continue)) {
            query = list(action = "query", prop = "revisions", format = "json", rvlimit = "max", rvprop = "user|timestamp|size", titles = page, 
                redirects = "")
        } else {
            query = list(action = "query", prop = "revisions", format = "json", rvlimit = "max", rvprop = "user|timestamp|size", titles = page, 
                rvcontinue = continue, redirects = "")
        }
        
        exec <- tryCatch(GET(paste("https://", domain, ".wikipedia.org/w/api.php", sep = ""), query = query), error = function(e) NULL)
        if (!is.null(exec)) {
            content <- content(exec, "parsed")
        } else {
            content <- NULL
        }
        continue <- tryCatch(content$continue$rvcontinue, error = function(e) NULL)
        
        content <- tryCatch(content[["query"]][[1]][[1]][[4]], error = function(e) NULL)
        
        if (!is.null(content)) {
            user <- sapply(content, function(x) {
                x$user
            })
            timestamp <- sapply(content, function(x) {
                x$timestamp
            })
            size <- sapply(content, function(x) {
                x$size
            })
            anon <- sapply(content, function(x) {
                if (is.null(x$anon)) 
                  FALSE else TRUE
            })
            data <- data.frame(user, timestamp, size, anon)
            if (length(data$size) > 1) {
                data$weight <- c(data[1:(nrow(data) - 1), 3] - data[2:nrow(data), 3], data[nrow(data), 3])
            } else {
                data$weight <- data$size
            }
            result <- rbind(result, data)
        }
        if (is.null(continue)) {
            break
        }
        
    }
    result <- result[-1, ]
    if (contrib) {
        result <- unique(result$user)
    }
    return(result)
}

#' Extracting the first contributors freom a list of revisions
#'
#' @param revisions A list of revisions built with page_revisions
#' @param threesold An integer, he function will return the \code{threesold} percent first contributors
#' @param keepAnon A boolean indicate weither the anonymous should be kept into the list or not
#'
#' @importFrom plyr ddply
#'
#' @return A character vector
#' @export
#' 
#' @family page functions
#'
#' @examples
#' # Downloading the list of revisions
#' revisions <- page_revisions('Action')
#' 
#' # Extracting the list of the 10 percent first contributor of the list
#' page_ranking(revisions,10,TRUE) 
page_ranking <- function(revisions, threesold, keepAnon = TRUE) {
    
    if (!keepAnon) {
        revisions <- as.data.frame(revisions[as.data.frame(revisions[, 4]) == FALSE, ])
    }
    
    revisions <- as.data.frame(revisions)
    data <- tryCatch(data.frame(name = revisions[, 1], weight = abs(as.numeric(revisions[, 5]))), error = function(e) NULL)
    if (!is.null(data)) {
        total <- sum(data[, 2])
        
        
        data <- ddply(data, ~name, function(x) {
            sum(x[, 2])
        })
        
        
        numberContrib <- as.integer((threesold * length(data[, 2]))/100)
        data <- data[order(data[, 2], decreasing = TRUE), ]
        result <- list(contrib = data[1:numberContrib, 1], weight = sum((data[1:numberContrib, 2]/total) * 100))
    } else {
        result <- FALSE
    }
    
    return(result)
}

#' Getting some quantitative infos about a wiki page
#'
#' @param title The title of the page
#' @param domain The domain where the wiki is located
#'
#' @return A named list
#' @export
#' 
#' @family page functions
#'
#' @examples
#' # Getting a list of information for the 'Nationalism' page in the french wiki
#' page_infos('Nationalisme') 
page_infos <- function(title, domain = "fr") {
    
    revisions <- page_revisions(title, domain)
    
    links <- page_links(title, domain = domain)
    
    ranking <- page_ranking(revisions, 10)
    
    return(list(`number of contributor` = length(unique(revisions$user)), `percent of anonymous contributions` = page_anon(revisions), `number of links from the article` = links$nb, 
        `number of revisions` = length(revisions[, 1]), `size in bytes` = sum(revisions$size), `percent of contribution due to the 10% first contributor` = ranking$weight))
    
} 
