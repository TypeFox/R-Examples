#'@title Retrieve Pageview Data for an Article
#'@description retrieves the pageview data for a particular article on a project, within
#'a provided time-range.
#'
#'@param project the name of the project, structured as \code{[language_code].[project]}
#'(see the default).
#'
#'@param article the article(s) you want to retrieve data for. Ideally features underscores in the title
#'instead of spaces, but happily converts if you forget to do this.
#'
#'@param platform The platform the pageviews came from; one of "all", "desktop", "mobile-web" and
#'"mobile-app". Set to "all" by default.
#'
#'@param user_type the type of users. One of "all", "user", "spider" or "bot". "all" by default.
#'
#'@param start the start \code{YYYYMMDDHH} of the range you want to cover. This can be
#'easily grabbed from R date/time objects using \code{\link{pageview_timestamps}}.
#'
#'@param end the end \code{YYYYMMDDHH} of the range you want to cover. NULL by default, meaning
#'that it returns 1 day of data.
#'
#'@param reformat Whether to reformat the results as a \code{\link{data.frame}} or not. TRUE by default.
#'
#'@param ... further arguments to pass to httr's GET.
#'
#'@seealso \code{\link{top_articles}} for the top articles per project in a given date range,
#'and \code{\link{project_pageviews}} for per-project pageviews.
#'
#'@examples
#'# Basic example
#'r_pageviews <- article_pageviews()
#'
#'# Modify the article
#'obama_pageviews <- article_pageviews(article = "Barack_Obama")
#'
#'@export
article_pageviews <- function(project = "en.wikipedia", article = "R (programming language)",
                              platform = "all", user_type = "all",
                              start = "2015100100", end = NULL, reformat = TRUE, ...){
  
  article <- gsub(x = article, pattern = " ", replacement = "_", fixed = TRUE)
  if(length(article) > 1){
    return(lapply(article, article_pageviews, project = project, platform = platform,
                  user_type = user_type, start = start, end = end, ...))
  }
  # Handle timestamps
  if(is.null(end)){
    end <- start
  }
  
  # Construct parameters
  parameters <- paste("per-article", project, ifelse(platform == "all", "all-access", platform),
                      ifelse(user_type == "all", "all-agents", user_type), article, "daily",
                      start, end, sep = "/")
  
  # Run and return
  data <- pv_query(parameters, ...)$items
  
  if(reformat){
    nameset <- names(data[[1]])
    data <- data.frame(matrix(unlist(data), nrow = length(data), byrow = TRUE), stringsAsFactors = FALSE)
    names(data) <- nameset
    data$views <- as.numeric(data$views)
    data <- data[,!names(data) == "granularity"]
  }
  return(data)
}

#'@title Retrieve Data on Top Articles
#'@description \code{top_articles} grabs data on the top articles for a project
#'in a given time period, and for a particular platform.
#'
#'@param project the name of the project, structured as \code{[language_code].[project]}
#'(see the default).
#'
#'@param platform The platform the pageviews came from; one of "all", "desktop", "mobile-web" and
#'"mobile-app". Set to "all" by default.
#'
#'@param year The year the articles were "top" in. 2015 by default.
#'
#'@param month The month the articles were "top" in. "10" by default; can be set to "all", for all
#'the months in \code{year}. If so, \code{day} must also be "all".
#'
#'@param day The day the articles were "top" in. "01" by default; can be set to "all", for all
#'the days in \code{month}.
#'
#'@param reformat Whether to reformat the results as a \code{\link{data.frame}} or not. TRUE by default.
#'
#'@param ... further arguments to pass to httr's GET.
#'
#'@seealso \code{\link{article_pageviews}} for per-article pageviews and \code{\link{project_pageviews}} for
#'per-project pageviews.
#'
#'@examples
#'# Basic example
#'enwiki_top_articles <- top_articles()
#'
#'# Use a narrower platform
#'enwiki_mobile_top <- top_articles(platform = "mobile-web")
#'
#'@importFrom jsonlite fromJSON
#'@export
top_articles <- function(project = "en.wikipedia", platform = "all", year = "2015",
                         month = "10", day = "01", reformat = TRUE, ...) {

  parameters <- paste("top", project, ifelse(platform == "all", "all-access", platform),
                      year, ifelse(month == "all", "all-months", month),
                      sep = "/", ifelse(day == "all", "all-days", day))
  results <- pv_query(parameters, ...)$items
  
  if(reformat){
    results <- do.call("rbind", lapply(results, function(x){
      
      # Reformat the list as a data.frame
      data <- x$articles
      reserved_names <- names(data[[1]])
      data <- data.frame(matrix(unlist(data), nrow = length(data), byrow = TRUE), stringsAsFactors = FALSE)
      names(data) <- reserved_names
      
      # Retype
      data$views <- as.numeric(data$views)
      data$rank <- as.numeric(data$rank)
      
      # Add new fields
      data$project <- x$project
      data$access <- x$access
      data$year <- as.integer(x$year)
      data$month <- as.integer(x$month)
      data$day <- as.integer(x$day)
      
      # Return
      return(data)
    }))
  }
  return(results)
}

#'@title Retrieve Per-Project Pageview Counts
#'
#'@description Retrieve pageview counts for a particular project.
#'
#'@param project the name of the project, structured as \code{[language_code].[project]}
#'(see the default).
#'
#'@param platform The platform the pageviews came from; one of "all", "desktop", "mobile-web" and
#'"mobile-app". Set to "all" by default.
#'
#'@param user_type the type of users. One of "all", "user", "spider" or "bot". "all" by default.
#'
#'@param granularity the granularity of data to return; do you want hourly or daily counts? Set
#'to "daily" by default.
#'
#'@param start the start \code{YYYYMMDDHH} of the range you want to cover. This can be
#'easily grabbed from R date/time objects using \code{\link{pageview_timestamps}}
#'
#'@param end the end \code{YYYYMMDDHH} of the range you want to cover. NULL by default, meaning
#'that it returns 1 day/hour of data (depending on the value passed to \code{granularity}).
#'
#'@param reformat Whether to reformat the results as a \code{\link{data.frame}} or not. TRUE by default.
#'
#'@param ... further arguments to pass to httr's GET.
#'
#'@examples
#'# Basic call
#'enwiki_1_october_pageviews <- project_pageviews()
#'
#'# Break it down to hourly
#'enwiki_hourly <- project_pageviews(granularity = "hourly", end = "2015100123")
#'
#'@seealso \code{\link{top_articles}} for the top articles per project in a given date range,
#'and \code{\link{article_pageviews}} for per-article pageviews.
#'
#'@export
project_pageviews <- function(project = "en.wikipedia", platform = "all", user_type = "all",
                              granularity = "daily", start = "2015100100", end = NULL, reformat = TRUE,
                              ...){
  
  # Handle timestamps
  if(is.null(end)){
    end <- start
  }
  
  # Construct parameters
  parameters <- paste("aggregate", project, ifelse(platform == "all", "all-access", platform),
                      ifelse(user_type == "all", "all-agents", user_type), granularity,
                      start, end, sep = "/")
  
  # Run
  data <- pv_query(parameters, ...)$items
  
  # Reformat if necessary, return either way.
  if(reformat){
    nameset <- names(data[[1]])
    data <- data.frame(matrix(unlist(data), nrow = length(data), byrow = TRUE), stringsAsFactors = FALSE)
    names(data) <- nameset
    data$views <- as.numeric(data$views)
  }
  return(data)
}
