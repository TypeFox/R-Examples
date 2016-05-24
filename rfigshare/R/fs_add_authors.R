#' Add author to an article
#' 
#' @author Carl Boettiger \email{cboettig@@gmail.com}
#' @param article_id id number of an article on figshare 
#' @param  authors a list or character string of authors or author id numbers (or mixed).  
#' @param session (optional) the authentication credentials from \code{\link{fs_auth}}. 
#' If not provided, will attempt to load from cache as long as figshare_auth has been run.  
#' @param create_missing (logical) Attempt to create authors not already registered on FigShare? 
#' @param debug return the httr result visibly?  
#' @return adds the requested authors to the given article
#' @export
#' @examples \dontrun{
#'  fs_add_authors("138", list("Scott Chamberlain", "Karthik Ram"))
#'  fs_add_authors("138", c("Scott Chamberlain", "Karthik Ram"))
#'  fs_add_authors("138", list("Scott Chamberlain", "97306"))
#'  fs_add_authors("138", list("Scott Chamberlain", 97306))
#'  fs_add_authors(138, 97306)
#' } 
fs_add_authors  <- function(article_id, authors, 
                            session = fs_get_auth(), 
                            create_missing = TRUE,
                            debug = FALSE){

  # See if any authors have been given as numeric id numbers instead of names:
  already_numeric <- sapply(authors, function(author) suppressWarnings(!is.na(as.numeric(author))) )
  # Go ahead and add those authors
  results <- sapply(authors[already_numeric], function(author_id) fs_add_author(article_id, author_id, session))
  # Get the IDs of the authors given
  remaining_authors <- authors[!already_numeric]
  if(length(remaining_authors)>0){ # don't try if all authors were given as ID numbers already 
    ids <- fs_author_ids(remaining_authors,  session)
    absent <- sapply(ids, is.null)
    missing_authors <- authors[absent]
    # Register an ID for any missing authors
    if(create_missing){
      missing_ids <- sapply(missing_authors, 
           function(author) fs_create_author(author, session))
      ids[absent] <- missing_ids
    } else if(any(absent)) {
      warning(paste0("Cannot find accounts for ", 
                     missing_authors, 
                     " so they cannot be added", collapse="\n"))
      ids <- ids[!absent]
    } else {
      message("found ids for all authors")
    }
    
    if (length(ids) == 0){
      warning(paste0("No matches found for authors", authors, collapse="\n"))
    }
    ## add the authors by ID number
    additional_results <- sapply(ids, function(author_id) fs_add_author(article_id, author_id, session))
    results <- c(results, additional_results)
    invisible(ids)
  } else {
    invisible(authors)
  }
  if(debug)
    results
  else
    invisible(results)
}


#' Get Author IDs from names
#' 
#' Take an author list, search for each author and return their FigShare ID.  
#' If no author is found, call fs_create_author and return the ID.
#' If multiple matches are found, allow user to choose interactively
#' @param authors a list/vector of authors (not a character string)
#' @param graphics logical (default False) use graphic input to disambiguate?
#' @param session (optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache. 
#' @return a list of author id numbers, or NULLS where ids cannot be found.  
#' @keywords internal
fs_author_ids <- function(authors, session = fs_get_auth(), graphics=FALSE){
  authors_info <- lapply(authors, function(author){
                         matches <- fs_author_search(author, session)
                         if(length(matches) == 0){
                           out <- list("author not found") 
                         } else if(length(matches) > 1){
                           opts <- c(sapply(matches, `[[`, "full_name"), "None of these")
                           a <- select.list(opts, 
                                  title = paste("Multiple matches found for", author, 
                                         "select which one you want"),
                                            graphics=graphics)
                           i <- which(opts == a)
                           if(i>length(matches))
                              out <- list("author not found")
                           else 
                             out <- matches[[i]]
                         } else if(length(matches) == 1){
                           out <- matches[[1]] 
                        }
                        out 
                  })
  ids <- sapply(authors_info, `[[`, "id")
}



#' Add an author to an article by ID number
#'
#' @param article_id id number of an article on figshare 
#' @param author_id the id number of a registered figshare user (see \code{\link{fs_author_search}})
#' @param session (optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache.  
#' @import RJSONIO
#' @keywords internal
fs_add_author <- function(article_id, author_id, session = fs_get_auth()){
  base <- "http://api.figshare.com/v1"
  method <- paste("my_data/articles", article_id, "authors", sep= "/")
  request = paste(base, method, sep="/")
  body <- toJSON(list("author_id"=author_id))
  config <- c(config(token = session), 
              add_headers("Content-Type" = "application/json"))
  out <- PUT(request, config=config,  body=body)
}

