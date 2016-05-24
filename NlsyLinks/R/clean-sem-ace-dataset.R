#' @name CleanSemAceDataset
#' @export
#' 
#' @title Produces a cleaned dataset that works well with when using SEM to estimate a univariate ACE model.
#' 
#' @description This function takes a `GroupSummary' \code{data.frame} (which is created by the \code{RGroupSummary} function) and returns a \code{data.frame} that is used by the \code{Ace} function.
#' @usage CleanSemAceDataset(dsDirty, dsGroupSummary, oName_S1, oName_S2, rName = "R")
#' 
#' @param dsDirty This is the \code{data.frame} to be cleaned.
#' @param dsGroupSummary The \code{data.frame} containing information about which groups should be included in the analyses.  It should be created by the \code{RGroupSummary} function.
#' @param oName_S1 The name of the manifest variable (in \code{dsDirty}) for the first subject in each pair.
#' @param oName_S2 The name of the manifest variable (in \code{dsDirty}) for the second subject in each pair.
#' @param rName The name of the variable (in \code{dsDirty}) indicating the pair's relatedness coefficient.
#'
#' @details The function takes \code{dsDirty} and produces a new \code{data.frame} with the following features:
#' 
#' [A] Only three existing columns are retained: \code{O1}, \code{O2}, and \code{R}.  They are assigned these names.
#' 
#' [B] A new column called \code{GroupID} is created to reflect their group membership (which is based on the \code{R} value).  These valuesa re sequential integers, starting at 1.  The group with the weakest \code{R} is 1.  The group with the strongest \code{R} has the largest \code{GroupID} (this is typically the MZ tiwns).
#' 
#' [C] Any row is excluded if it has a missing data point for \code{O1}, \code{O2}, or \code{R}.
#' 
#' [D] The \code{data.frame} is sorted by the \code{R} value.  This helps program against the multiple-group SEM API sometimes.
#'
#'@return A \code{data.frame} with one row per subject pair.  The \code{data.frame} contains the following variables (which can NOT be changed by the user through optional parameters):
#' \item{R}{The pair's \code{R} value.}
#' \item{O1}{The outcome variable for the first subject in each pair.}
#' \item{O2}{The outcome variable for the second subject in each pair.}
#' \item{GroupID}{ Indicates the pair's group membership.}
#'
#' @author Will Beasley
#' 
#' @examples
#' library(NlsyLinks) #Load the package into the current R session.
#' dsLinks <- Links79PairExpanded #Start with the built-in data.frame in NlsyLinks
#' dsLinks <- dsLinks[dsLinks$RelationshipPath=='Gen2Siblings', ] #Use only NLSY79-C siblings
#' 
#' oName_S1 <- "MathStandardized_S1" #Stands for Outcome1
#' oName_S2 <- "MathStandardized_S2" #Stands for Outcome2
#' dsGroupSummary <- RGroupSummary(dsLinks, oName_S1, oName_S2)
#' 
#' dsClean <- CleanSemAceDataset( dsDirty=dsLinks, dsGroupSummary, oName_S1, oName_S2, rName="R" )
#' summary(dsClean)
#' 
#' dsClean$AbsDifference <- abs(dsClean$O1 - dsClean$O2)
#' plot(jitter(dsClean$R), dsClean$AbsDifference, col="gray70")
#' @keywords ACE


CleanSemAceDataset <- function( dsDirty, dsGroupSummary, oName_S1, oName_S2, rName="R" ) {
  rLevelsToInclude <- dsGroupSummary[dsGroupSummary$Included, rName]
  
  #It's necessary to drop the missing Groups & unnecessary columns.  Missing O1s & O2s are dropped for the sake of memory space.
  oldColumnNames <- c(rName, oName_S1, oName_S2)
  newColumnNames <- c("R", "O1", "O2")
  selectedRows <- (
    (!base::is.na(dsDirty[, rName])) & 
    (dsDirty[, rName] %in% rLevelsToInclude) & 
    (!base::is.na(dsDirty[, oName_S1])) & 
    (!base::is.na(dsDirty[, oName_S2]))
  )
  
  dsClean <- dsDirty[selectedRows, oldColumnNames] 
  
  colnames(dsClean) <- newColumnNames
  
  dsClean <- dsClean[base::order(dsClean$R), ] #TODO: Rewrite overall code so this statement is not longer necessary anyomre.
  
  #This helper function allows for slight imprecision from floating-point arithmetic.
  EqualApprox <- function( target, current, toleranceAbsolute=1e-8) {  
    return( abs(target-current) < toleranceAbsolute ) 
  }
  
  #rLevelsToExclude <- dsGroupSummary[!dsGroupSummary$Included, 'R']
  
  #This loop assigns a GroupID, depending on their R value. TODO: possibly rewrite and vectorize with plyr.
  dsClean$GroupID <- NA
  for( groupIndex in seq_along(rLevelsToInclude) ) {
    r <- rLevelsToInclude[groupIndex]
    memberIndices <- base::sapply(dsClean$R, EqualApprox, r)
    dsClean$GroupID[memberIndices] <- groupIndex
  }
  
  return( dsClean )
}
