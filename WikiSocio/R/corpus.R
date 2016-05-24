#' Creating a page corpus
#'
#' @param method The method employed to create the corpus. \describe{
#'  \item{category}{Every page in the category x is included to the corpus}
#'  \item{random}{Select x random page}
#'  \item{xpath}{x is a vector containing in first position the URL of one page,
#'  and in second position the wpath request to apply to this URL.
#'  Will return a list of words, wich are the names of pages if the xpath request is writed correctly}
#' }
#' @param x The pointer of the method.
#' @param limit If not 'max', then an integer giving the length of the corpus. Useless for random method.
#' @param domain The domain where the wiki is located
#'
#' @importFrom httr GET content
#' @importFrom stringr str_split_fixed
#'
#' @family corpus functions
#'
#' @return A character vector.
#' @export
#'
#' @seealso \code{\link{corpus_contrib_create}} for contributors corpus, \code{\link{corpus_page_data}} to put data into the corpus you made with this function
#'
#' @examples
#' # Creating a page corpus formed by 3 random selected page
#' corpus_page_create('random',3)
corpus_page_create <- function(method, x, limit = "max", domain = "fr") {
    methods <- c("category", "random", "xpath")
    if (!is.element(method, methods)) 
        stop("Please enter a valid method name")
    
    if (method == "category") {
        continue <- NULL
        corpus <- matrix(ncol = 2)
        repeat {
            if (is.null(continue)) {
                query = list(action = "query", list = "categorymembers", format = "json", cmtitle = x, cmprop = "title", cmlimit = "max")
            } else {
                query = list(action = "query", list = "categorymembers", format = "json", cmtitle = x, cmprop = "title", cmlimit = "max", cmcontinue = continue)
            }
            
            
            exec <- GET(paste("https://", domain, ".wikipedia.org/w/api.php", sep = ""), query = query)
            content <- content(exec, "parsed")
            continue <- content$continue$cmcontinue
            
            data <- tryCatch(content$query$categorymembers, error = function(e) NULL)
            if (!is.null(data)) {
                data <- matrix(unlist(data), ncol = 2, byrow = TRUE)
                corpus <- rbind(corpus, data)
            }
            if (is.null(continue)) {
                break
            }
            
        }
        
        corpus <- corpus[-1, 2]
    }
    if (method == "random") {
        corpus <- random(x, domain = domain)
    }
    if (method == "xpath") {
        corpus <- xpath(x[1], x[2])
        if (is.integer(limit)) {
            if (limit < length(corpus)) {
                corpus <- corpus[1:limit]
            }
        }
    }
    return(corpus)
}

#' Adding page data to a page corpus
#'
#' @param method The method employed to get data into the corpus. \describe{
#'  \item{variables}{For each page of the corpus, return variables specified in selection}
#' }
#' @param x A character vector created with \code{\link{corpus_page_create}}
#' @param selection A character vector giving all the data variables to get in the corpus : \describe{
#'  \item{nbLinks}{Number of links in the page}
#'  \item{nbContrib}{Number of contributors who have edited the page at least once}
#'  \item{nbRevisions}{Number of revisions of the page}
#'  \item{percentAnon}{Percentage of anonymous contributions}
#'  \item{percent10}{Percentage of the text due to the 10 \% first contributor}
#' }
#' @param domain The domain where the wiki is located
#'
#' @import pbapply
#'
#' @return A data-frame.
#' @export
#'
#' @family corpus functions
#'
#'
#' @examples
#' # Creating a page corpus with 3 randomly selected page
#' corpus <- corpus_page_create('random',3)
#' corpus <- corpus_page_data("variables",corpus)
corpus_page_data <- function(method, x, selection = c("nbLinks", "nbContrib", "nbRevisions", "percentAnon", "percent10"), domain = "fr") {
    methods <- c("variables")
    if (!is.element(method, methods)) 
        stop("Please enter a valid method name")
    
    if (method == "variables") {
        corpus <- data.frame(article = x)
        if ("nbLinks" %in% selection) {
            print("Number of links")
            corpus$nbLiens <- pbsapply(as.character(corpus$article), function(x) {
                return(tryCatch(page_links(x, domain = domain)$nb, error = function(e) NULL))
            })
        }
        
        print("Revisions' list extraction")
        list <- pblapply(corpus$article, function(x) {
            return(tryCatch(page_revisions(as.character(x), domain = domain), error = function(e) NULL))
        })
        names(list) <- corpus$article
        
        if ("nbContrib" %in% selection) {
            print("Number of contributors")
            corpus$nbContrib <- pbsapply(as.character(corpus$article), function(x) {
                result <- unique(list[x][[1]]$user)
                return(tryCatch(length(result), error = function(e) NULL))
            })
        }
        
        if ("nbRevisions" %in% selection) {
            print("Number of revisions")
            corpus$nbRevisions <- pbsapply(corpus$article, function(x, list) {
                return(tryCatch(length(as.data.frame(list[as.character(x)])[, 1]), error = function(e) NULL))
            }, list)
        }
        
        if ("percentAnon" %in% selection) {
            print("Anonymous contributions percentage")
            corpus$percentAnon <- pbsapply(corpus$article, function(x, list) {
                data <- page_anon(list[as.character(x)])
                return(tryCatch(data, error = function(e) NULL))
            }, list)
        }
        
        if ("percent10" %in% selection) {
            print("Weigth of ten percent principal contributors")
            corpus$weightContrib <- pbsapply(corpus$article, function(x, list) {
                data <- page_ranking(list[as.character(x)], 10)
                if (is.null(data$weight)) {
                  return(NULL)
                } else {
                  return(data$weight)
                }
            }, list)
        }
        
    }
    
    return(corpus)
}

#' Creating a contributor corpus
#'
#' @param method The method employed to get data into the corpus. \describe{
#'   \item{category}{Every contributor page included in the category x would be
#'   included in the corpus} \item{random}{Select x random contributor page}
#'   \item{xpath}{x is a vector containing in first position the URL of one
#'   page, and in second position the wpath request to apply to this URL. Will
#'   return a list of words, wich are the names of the contributors if the xpath
#'   request is writed correctly} }
#' @param x The pointer of the method.
#' @param limit If not 'max', then an integer giving the corpus size. Useless
#'   for random method.
#' @param domain The domain where the wiki is located.
#'
#' @importFrom httr GET content
#' @importFrom stringr str_split_fixed
#' @import pbapply
#'
#' @return A character vector
#' @export
#'
#' @family corpus functions
#'
#' @seealso \code{\link{corpus_contrib_data}} to put data into the corpus you
#'   made with this function, \code{\link{corpus_contrib_select}} to clean the
#'   corpus with a specific criteria
#'
#' @examples
#' # Creating a contributor corpus formed by contributors of the 'Action' page on the french wiki.
#' corpus_contrib_create('page','Action')
#'
#' # Same as previous, by limiting the size of the corpus at 5 contributors.
#' corpus_contrib_create('page','Action',limit=5)
corpus_contrib_create <- function(method, x, limit = "max", domain = "fr") {
    methods <- c("category", "random", "xpath", "page")
    
    if (!is.element(method, methods)) 
        stop("Please enter a valid method name")
    
    if (method == "category") {
        corpus <- as.data.frame(corpus_page_create("category", x, limit, domain = domain))
        corpus <- corpus[corpus[, 1] == "2", ]
        
        corpus <- str_split_fixed(corpus[, 2], ":", 2)[, 2]
        corpus <- str_split_fixed(corpus, "/", 2)[, 1]
    }
    
    if (method == "page") {
        corpus <- page_revisions(x, contrib = TRUE, domain = domain)
    }
    
    if (method == "random") {
        corpus <- random(x, namespace = "2", domain = domain)
        corpus <- str_split_fixed(corpus, "/", 2)[, 1]
        corpus <- str_split_fixed(corpus, ":", 2)[, 2]
    }
    
    if (method == "xpath") {
        corpus <- xpath(x[1], x[2])
        if (is.integer(limit)) {
            if (limit < length(corpus)) {
                corpus <- corpus[1:limit]
            }
        }
    }
    
    return(corpus)
    
}

#' Adding data to a contributor corpus
#'
#' @param method The method employed to get data into the corpus. \describe{
#'  \item{raw}{For each user, downloading caracteristics and testing other member of the corpus.}
#'  \item{page}{For each user, downloading edited pages.}
#' }
#' @param x The users' corpus.
#' @param prefix In the wiki, the prefix for users' pages namespace.
#' @param domain The domain where the wiki is located.
#'
#' @import pbapply
#' @importFrom stringr str_replace
#'
#' @family corpus functions
#'
#' @return A data-frame.
#' @export
#'
#' @examples
#' c <- corpus_contrib_create('page','Action',limit=5)
#' \donttest{
#'  c <- corpus_contrib_data('raw',c)
#' }
corpus_contrib_data <- function(method, x, prefix = "Utilisateur", domain = "fr") {
    methods <- c("raw", "page")
    if (!is.element(method, methods)) 
        stop("Please enter a valid method name")
    
    if (method == "raw") {
        print("Downloading user caracteristics..")
        corpus <- pblapply(unique(x), function(x) {
            cat <- page_cat(paste(prefix, x, sep = ":"), domain = domain)
        })
        
        names(corpus) <- unique(x)
        
        dimDF <- unique(unlist(corpus))
        
        print("Formatting data table...")
        corpus <- pblapply(corpus, function(x, dimDF) {
            as.vector(sapply(dimDF, function(y, x) {
                if (y %in% x) 
                  TRUE else FALSE
            }, x))
        }, dimDF)
        
        corpus <- as.data.frame(do.call(rbind, corpus))
        
        corpus <- list(index = dimDF, data = corpus)
    }
    
    if (method == "page") {
        print("Extracting the list of contributions done by users into corpus")
        
        corpus <- pblapply(x, function(x) {
            page <- contrib_list(x, page = TRUE, domain = domain)
            data.frame(contrib = rep(x, length(page)), article = page)
        })
        
        corpus <- do.call(rbind, corpus)
    }
    
    return(corpus)
}

#' Selecting the member of a corpus, based on their contributions.
#'
#' This function is used to have a corpus of big contributors of the wiki.
#' To use this, you need a contributor corpus completed with 'page' method of corpus_contrib_data
#'
#' BE CAUTIOUS : this function is VERY time-consuming.
#'
#' @param method Method employed to get data into the corpus. \describe{
#'  \item{firstContrib}{Select only contributors who are part of the \code{threesold} \% of first contibutor of each page}
#' }
#' @param x A corpus created with corpus_contrib_create
#' @param threesold An integer used as a threesold to decide weither a contributor can be on the corpus or not. For instance, if \code{threesold=5}, only contribor who are is in the 5 percent first contributor of the article.
#' @param domain The domain where the wiki is located
#'
#' @import pbapply
#'
#' @return A data-frame.
#' @export
#'
#' @family corpus functions
#'
#' @examples
#' # creating a corpus of 5 contributor of the action page
#' c <- corpus_contrib_create('page','Action',limit=5)
#' \donttest{
#' c <- corpus_contrib_data('page',c)
#'
#' # Keeping on this corpus only the contributor who are part of
#' # the 5 percent first contributors of the correspondant article.
#' c <- corpus_contrib_select(c,5)
#' }

corpus_contrib_select <- function(method, x, threesold, domain = "fr") {
    methods <- c("firstContrib")
    if (!is.element(method, methods)) 
        stop("Please enter a valid method name")
    
    if ((!is.data.frame(x)) && names(x) == c("contrib", "article")) 
        stop("Please enter a corpus made with corpus_contrib_data with the method page - after creating him with corpus_contrib_create.")
    
    if (method == "firstContrib") {
        corpus <- x
        
        print("Revisions' list extraction...")
        
        list <- pblapply(unique(corpus$article), function(x) {
            page_revisions(as.character(x), domain = domain)
        })
        names(list) <- unique(corpus$article)
        
        print("First contributors' list extraction...")
        
        rank <- pblapply(unique(corpus$article), function(x, list, threesold) {
            page_ranking(list[as.character(x)], threesold)
        }, list, threesold)
        names(rank) <- unique(corpus$article)
        
        print("Corpus cleaning...")
        
        corpus$weight <- pbapply(corpus, 1, function(x, rank) {
            x[1] %in% rank[[as.character(x[2])]]$contrib
        }, rank)
        
        corpus <- corpus[corpus$weight == TRUE, c(1,2)]
        
    }
    
    return(corpus)
    
}


#' @importFrom httr GET content

random <- function(nb, domain = "fr", namespace = "0") {
    quotient <- nb%/%500
    reste <- nb%%500
    
    result <- vector()
    
    if (quotient > 0) {
        for (i in 1:quotient) {
            query <- list(action = "query", list = "random", format = "json", rnlimit = "max", rnnamespace = namespace, redirects = "")
            
            exec <- GET(paste("https://", domain, ".wikipedia.org/w/api.php", sep = ""), query = query)
            content <- content(exec, "parsed")
            content <- content[[4]][[1]]
            content <- sapply(content, function(x) {
                x$title
            })
            result <- c(result, content)
        }
    }
    
    if (reste != 0) {
        query <- list(action = "query", list = "random", format = "json", rnlimit = reste, rnnamespace = namespace, redirects = "")
        
        exec <- GET(paste("https://", domain, ".wikipedia.org/w/api.php", sep = ""), query = query)
        content <- content(exec, "parsed")
        content <- content[["query"]][[1]]
        content <- sapply(content, function(x) {
            x$title
        })
        
        result <- c(result, content)
        result <- unique(result)
        
        if (length(result) != nb) {
            print(paste("Vector is shorter than ", nb, " . ", length(result), " pages really selected.", sep = ""))
        }
    }
    return(result)
} 
