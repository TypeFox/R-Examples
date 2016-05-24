#' Assigns title/abstract screening efforts to a team.  
#'
#' Randomly distributes screening tasks evenly or unevenly across multiple team
#' members.  It populates this effort in a data frame column that includes this
#' screening work (e.g., ABSTRACTS and TITLES).    
#'
#' @param aDataFrame A data.frame containing the titles and abstracts to be 
#'    screened by a team.  The default assumes that the data.frame has already
#'    been formatted using \code{effort_initialize}.  This data.frame will be
#'    populated with screening efforts.  See example: 
#'    \code{\link{example_references_metagear}}
#' @param dual When \code{TRUE}, distributes effort using a dual screening 
#'    design where two members will screen the same random collection of 
#'    titles/abstracts.  Requires the team to have an even number of members.
#' @param reviewers A vector with the names of each team member.
#' @param column_name Changes the default label of the "REVIEWERS" column 
#'    that contains the screening efforts of each team member.
#' @param effort A vector of percentages used to allocate screening
#'    tasks among each team member.  When not called explicitly, assumes effort 
#'    to be distributed evenly among all members.  Must be the same length as 
#'    the number of team members, and also sum to 100. 
#' @param initialize When \code{TRUE}, initializes the data.frame so that 
#'   efforts could be distributed, calls: \code{\link{effort_initialize}}.  
#'    Default is \code{FALSE}.
#' @param save_split Saves the allocated team effort into separate effort_*.csv
#'    files for individual screening tasks.  These files can be given to each
#'    member to screen their random title/abstract subset.  All files can be 
#'    merged once all screening tasks have been completed using 
#'    \code{\link{effort_merge}}.  
#' @param directory Changes the default location/directory for where the  
#'    effort_*.csv will be saved.  If not explicitly called, it will deposit 
#'    files in the current working directory.
#'
#' @return A data.frame with title/abstract screening efforts randomly 
#'    distributed across a team.
#'
#' @examples
#' data(example_references_metagear)
#' theTeam <- c("Christina", "Luc")
#' effort_distribute(example_references_metagear, initialize = TRUE, reviewers = theTeam)
#'
#' @seealso \code{\link{effort_initialize}}, \code{\link{effort_merge}}, 
#'    \code{\link{effort_summary}}
#'
#' @export effort_distribute

effort_distribute <- function (aDataFrame,
                               dual = FALSE,
                               reviewers, 
                               column_name = "REVIEWERS", 
                               effort = NULL, 
                               initialize = FALSE,
                               save_split = FALSE,
                               directory = getwd() ) {
  
  number_REFS <- nrow(aDataFrame)
  number_reviewers <- length(reviewers)
  
  # add REVIEWER, STUDY_ID, and INCLUDE columns to reference library
  if(initialize == TRUE) aDataFrame <- effort_initialize(aDataFrame,
                                                         dual = dual)
  
  if(dual == TRUE) {
    if(!is.null(effort)) .metagearPROBLEM("error", 
                                          "can only assign dual effort evenly among reviewers")
    if(number_reviewers %% 2 == 1) .metagearPROBLEM("error", 
                                                    "can only assign dual effort with an even number of reviewers")
    
    reviewers_A <- reviewers[1:number_reviewers %% 2 == 1]
    theEffort_A <- gl(length(reviewers_A), 
                      ceiling(number_REFS / length(reviewers_A)), 
                      number_REFS, 
                      labels = reviewers_A)

    reviewers_B <- reviewers[1:number_reviewers %% 2 == 0]
    theEffort_B <- gl(length(reviewers_B), 
                      ceiling(number_REFS / length(reviewers_B)), 
                      number_REFS, 
                      labels = reviewers_B)
    
    dualEffort <- data.frame(A = theEffort_A, B = theEffort_B)
    theEffort <- dualEffort[sample(nrow(dualEffort), 
                                   nrow(dualEffort), replace = FALSE), ]
    
    aDataFrame["REVIEWERS_A"] <- theEffort["A"]
    aDataFrame["REVIEWERS_B"] <- theEffort["B"]
  }
  else {
    
    # generate reviewers tasks evenly or via custom 'effort' 
    if(is.vector(effort) && length(unique(effort)) != 1) {
      if(sum(effort) != 100) .metagearPROBLEM("error", 
                                              "Effort does not sum to 100%.")
      theEffort <- rep(reviewers, round((number_REFS * (effort / 100))))
    } else {
      theEffort <- gl(number_REFS, 
                      ceiling(number_REFS / number_reviewers), 
                      number_REFS, 
                      labels = reviewers)
    }

    # randomly populate REVIEWERS column with tasks
    aDataFrame[, column_name] <- sample(theEffort, 
                                      length(theEffort), 
                                      replace = FALSE)
  }
  
  # splits reference library into seperate reviewer csv files and
  # hides teams if dual reviewing
  if(save_split == TRUE) {
    if(dual == TRUE) {
      removeVars <- names(aDataFrame) %in% c("REVIEWERS_B", "INCLUDE_B")
      effort_save(aDataFrame[!removeVars], 
                  column_name = "REVIEWERS_A", directory)
      removeVars <- names(aDataFrame) %in% c("REVIEWERS_A", "INCLUDE_A")
      effort_save(aDataFrame[!removeVars], 
                  column_name = "REVIEWERS_B", directory)
    }
    else effort_save(aDataFrame, column_name, directory)
  }
  
  return(aDataFrame)
}