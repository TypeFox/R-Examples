#' Formats a reference dataset for title/abstract screening efforts.
#'
#' Adds columns with standardized labels to a data framw with bibliographic data
#' on journal articles.  These columns will be used to assign reviewers, 
#' implementation of dual screening design, and the coding of
#' inclusion/exclusions screening decisions.    
#'
#' @param aDataFrame A data.frame object that includes the titles and 
#'    abstracts to be screened.  It will be formatted for screening efforts.  
#'    See example: \code{\link{example_references_metagear}}
#' @param study_ID When \code{FALSE}, does not add a column "STUDY_ID" that
#'    includes a unique identification number for each reference (row) in the
#'    dataFrame.
#' @param unscreenedValue Changes the default coding (a string) of "not vetted"
#'    that designates whether an abstract remains to be screened or vetted as
#'    part of the "INCLUDE" column.
#' @param dual When \code{TRUE}, formats dataFrame for a dual screening (paired)
#'    design.  Creates two reviewer teams: REVIEWERS_A and REVIEWERS_B.    
#' @param front When \code{FALSE}, adds new columns to the back end of the
#'    dataframe.  When \code{TRUE}, adds columns to the front.     
#'
#' @return A data.frame formatted for title/abstract screening efforts.
#'
#' @examples
#' data(example_references_metagear)
#' effort_initialize(example_references_metagear)        
#'
#' @seealso \code{\link{effort_distribute}}, \code{\link{effort_merge}}, 
#'    \code{\link{effort_summary}}
#'
#' @export effort_initialize

effort_initialize <- function(aDataFrame,
                              study_ID = TRUE,
                              unscreenedValue = "not vetted",
                              dual = FALSE, 
                              front = TRUE) {
  
  if(dual == TRUE) {
    new_cols <- data.frame(REVIEWERS_A = NA, 
                           INCLUDE_A = unscreenedValue,
                           REVIEWERS_B = NA, 
                           INCLUDE_B = unscreenedValue)
  }
  else new_cols <- data.frame(REVIEWERS = NA, INCLUDE = unscreenedValue)
  
  if(study_ID == TRUE) new_cols <- cbind(data.frame(STUDY_ID = 1:nrow(aDataFrame)), 
                                         new_cols)
  
  if(front == TRUE) newDataFrame <- cbind(new_cols, aDataFrame)
  else newDataFrame <- cbind(aDataFrame, new_cols)
  
  return(newDataFrame)
}
