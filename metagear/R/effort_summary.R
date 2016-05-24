#' Provides a text summary of screening efforts among the reviewing team.  
#'
#' Summarizes the number of studies screened, which were identified to be
#' included/excluded from the project, as well as those with conflicting
#' agreement on their inclusion/exclusion.  If a dual (paired) design was 
#' implemented to screen references, then it also provides inter-reviewer 
#' agreement estimate following Cohen (1960) that describes the agreement (or 
#' repeatability) of screening/coding decisions.     
#'
#' @param aDataFrame A data.frame containing the titles and abstracts that were 
#'    screened by a team.  The default assumes that the data.frame is the
#'    merged effort across the team using \code{effort_merge}.
#' @param column_reviewers Changes the default label of the "REVIEWERS" column 
#'    that contains the screening efforts of each team member.
#' @param column_effort Changes the default label of the "INCLUDE" column 
#'    that contains the screening decisions (coded references) of each team 
#'    member.
#' @param dual When \code{TRUE}, provides a summary of the dual screening
#'    effort as well as estimation of inter-reviewer agreements following 
#'    Cohen's (1960) kappa (K) and Landis and Koch's (1977) interpretation 
#'    benchmarks. 
#' @param quiet When \code{TRUE}, does not print to console the summary table. 
#'
#' @return A data frame with summary information on the screening tasks of a
#'    reviewing team.
#'
#' @seealso \code{\link{effort_initialize}}, \code{\link{effort_distribute}}, 
#'    \code{\link{effort_merge}} 
#'
#' @references Cohen, J. 1960. A coefficient of agreement for nominal scales. 
#'    Educational and Psychological Measurement 20: 37-46.
#' @references Landis, J.R., and Koch, G.G. 1977. The measurement of observer 
#'    agreement for categorical data. Biometrics 33: 159-174.
#'
#' @importFrom stats addmargins
#' @export effort_summary


effort_summary <- function(aDataFrame,
                           column_reviewers = "REVIEWERS",
                           column_effort = "INCLUDE",
                           dual = FALSE,
                           quiet = FALSE) {

  #TO DO: add error if data.frames have not beed screened
                           
  if(dual == TRUE) {
    aTable_A <- table(aDataFrame[, c("REVIEWERS_A", "INCLUDE_A")])
    aTable_A <- aTable_A[, order(unlist(dimnames(aTable_A)["INCLUDE_A"]))]
    aTable_A <- addmargins(aTable_A, FUN = list(TOTAL = sum), quiet = TRUE)
    names(dimnames(aTable_A)) <- c("", "")
    aTable_A <- as.data.frame.matrix(aTable_A)
    aTable_A["%"] <- (prop.table(aTable_A)*(nrow(aTable_A) * 100))[, nrow(aTable_A)]
    
    aTable_B <- table(aDataFrame[, c("REVIEWERS_B", "INCLUDE_B")])
    aTable_B <- aTable_B[, order(unlist(dimnames(aTable_B)["INCLUDE_B"]))]
    aTable_B <- addmargins(aTable_B, FUN = list(TOTAL = sum), quiet = TRUE)
    names(dimnames(aTable_B)) <- c("", "")
    aTable_B <- as.data.frame.matrix(aTable_B)
    aTable_B["%"] <- (prop.table(aTable_B)*(nrow(aTable_B) * 100))[, nrow(aTable_B)]
        
    theTeams <- split(aDataFrame, aDataFrame[, "REVIEWERS_A"])   
    theAggrements <- lapply(theTeams, function(x) {
                                        aTable <- table(x[, c("INCLUDE_A","INCLUDE_B")])
                                        aTable <- aTable[order(unlist(dimnames(aTable)["INCLUDE_A"])), order(unlist(dimnames(aTable)["INCLUDE_B"]))]
                                        return(as.matrix(as.data.frame.matrix(aTable)))
                                      }
                           )
    
    aggrees <- lapply(theAggrements, function(x) sum(diag(x)))
    total <- lapply(theAggrements, function(x) sum(x))
    rowSUMS <- lapply(theAggrements, function(x) rowSums(x))
    colSUMS <- lapply(theAggrements, function(x) colSums(x))
    
    for(i in 1:length(rowSUMS)) {
      ef <- sum((rowSUMS[[i]] * colSUMS[[i]]) / total[[i]])
      K <- (aggrees[[i]] - ef)/(total[[i]] - ef)
      if(i == 1) {
        theKappas <- K
        theScale <- Landis_Koch_scale(K)
      }
      else {
        theKappas <- c(theKappas, K)
        theScale <- c(theScale, Landis_Koch_scale(K))
      }
    }  
    
    theK <- paste0(c(theKappas, ""), c(rep(" (", length(theKappas)), ""), c(theScale, ""), c(rep(")", length(theKappas)), ""))
    
    aTable <- cbind("Team_A" = row.names(aTable_A), 
                    aTable_A, 
                    "   " = rep(" ", nrow(aTable_A)), 
                    "Team_B" = row.names(aTable_B), 
                    aTable_B, 
                    "   " = rep(" ", nrow(aTable_A)), 
                    "Cohen's K" = theK )
                    
    summarized <- Reduce('+', theAggrements)             
                    
  } else {
    aTable <- (addmargins(table(aDataFrame[, c(column_reviewers, column_effort)]), FUN = list(TOTAL = sum), quiet = TRUE))
    names(dimnames(aTable)) <- c("", "")
    aTable <- as.data.frame.matrix(aTable)
    aTable["%"] <- (aTable["TOTAL"] / aTable[nrow(aTable), "TOTAL"]) * 100
  }
  
  if(quiet == FALSE) {
    if(dual == TRUE) {
      cat("=== SCREENING EFFORT SUMMARY ===\n\n")
      cat(paste("  ", summarized[3,3], "candidate studies identified\n"))
      cat(paste("  ", summarized[2,2], "studies excluded\n"))
      cat(paste("  ", sum(summarized) - summarized[3,3] - summarized[2,2], "studies with conflicting agreement that need additional screening\n"))
      cat("  ----\n")
      cat(paste("  ", sum(summarized), "TOTAL SCREENED\n\n"))
      cat("=== DUAL SCREENING DESIGN SUMMARY ===\n\n")
      print(aTable, row.names = FALSE, right = FALSE)
      cat("\n")
    } else {
      cat("=== SCREENING EFFORT SUMMARY ===\n\n")
      cat(paste("  ", aTable["TOTAL", 3], "candidate studies identified\n"))
      cat(paste("  ", aTable["TOTAL", 2], "studies excluded\n"))
      cat(paste("  ", aTable["TOTAL", 1], "challenging studies needing additional screening\n"))
      cat("  ----\n")
      cat(paste("  ", aTable["TOTAL", 4], "TOTAL SCREENED\n\n"))
      cat("=== SCREENING DESIGN SUMMARY ===\n\n")
      print(aTable)
    }
  }
    
  return(aTable)
}