#' Get data from Statistics Netherlands (CBS)
#' 
#' Retrieves data from a table of Statistics Netherlands. A list of tables
#' can be retrieved with \code{\link{get_table_list}}. 
#' Optionaly the data can be filtered on category values. 
#' The filter is specified with \code{<column_name> = <values>} in which \code{<values>} is a character vector.
#' Rows with values that are not part of the character vector are not returned. Note that the values
#' have to be raw (un-recoded) values.
#' 
#' @note All data are downloaded using \code{\link{download_table}}
#' 
#' @param id Identifier of table, can be found in \code{\link{get_table_list}}
#' @param ... optional filter statemenets
#' @param recode Should the categories of the table be recoded with their title
#' (TRUE) or with their key (FALSE)? Should column names be coded with title (TRUE)
#' or key (FALSE) 
#' @param dir Directory where the table should be downloaded. Defaults to temporary
#' directory
#' @param base_url optionally specify a different server. Useful for
#' third party data services implementing the same protocol.
#' @return \code{data.frame} with the requested data. Note that a csv copy of the data is stored in \code{dir}.
#' @export
#' @examples 
#' \dontrun{
#' 
#' # get data for main (000000) Consumer Price Index (7196ENG) for March 2000, 
#'  get_data(id="7196ENG", Periods="2000MM03", CPI="000000")
#' }
get_data <- function(id, ..., recode=TRUE, dir=tempdir(), base_url = CBSOPENDATA){
  meta <- download_table(id, ..., dir=dir, cache=TRUE, base_url = base_url)
  data <- read.csv(file.path(dir, "data.csv"), colClasses="character")
  
  if (recode){
    dims <- names(data)[names(data) %in% names(meta)]
    for (d in dims){
      x <- as.factor(data[[d]])
      dim <- meta[[d]]
      levels(x) <- dim$Title[match(levels(x), dim$Key)]
      data[[d]] <- x
    }
    columns <- meta$DataProperties
    
    m <- match(columns$Key, colnames(data), nomatch = 0)
    
    # can cause duplicated names!
    #colnames(data)[m] <- columns$Title[m > 0]
    
    # TODO recode column names from meta$DataProperties and
    # convert columns to correct type.
  }
  
  class(data) <- c('tbl', 'tbl_df', 'data.frame')
  data
}

#get_data("81819NED")
