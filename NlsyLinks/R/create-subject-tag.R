#' @name CreateSubjectTag
#' @export
#' 
#' @title Creates a \code{SubjectTag}.  This value uniquely identifies subjects, when both generations are included in the same dataset.
#' @description A \code{SubjectTag} uniquely identify subjects.  For Gen2 subjects, the SubjectTag is identical to their CID (ie, C00001.00 -the SubjectID assigned in the NLSY79-Children files).  However for Gen1 subjects, the SubjectTag is their CaseID (ie, R00001.00), with "00" appended.  This manipulation is necessary to identify subjects uniquely in inter-generational datasets.  A Gen1 subject with an ID of 43 becomes 4300.  The SubjectTags of her four children remain 4301, 4302, 4303, and 4304.
#'
#' @usage CreateSubjectTag(subjectID, generation)
#' @param subjectID The ID assigned by the NLSY.  For Gen1 subjects, this will be their CaseID (ie, R00001.00).  For Gen2 subjects, this will be their CID (ie, C00001.00).
#' @param generation The generation of the subject.  Values are either 1 or 2, representing Gen1 and Gen2.  
#' 
#' @details For a fuller explanation of \code{SubjectTag} in context, see the \code{\link{Links79Pair}} dataset documentation.
#' 
#' @return A integer value under normal circumstances. An error is thrown if the vectors \code{subjectID} and \code{generation} are different lengths. If either input vector has NA values, the respective output element(s) will be NA too.
#' @author Will Beasley
#' @seealso \code{\link{Links79Pair}}
#' 
#' @examples
#' library(NlsyLinks) #Load the package into the current R session.
#'   
#' #Typically these two vectors will come from a data frame.
#' subjectIDs <- c(71:82, 10001:10012)
#' generation <- c(rep(1, 12), rep(2, 12))
#'   
#' CreateSubjectTag(subjectIDs, generation)
#' #Returns 7100, ..., 8200, 10001, ..., 10012
#'   
#' #Use the ExtraOutcomes79 dataset, with numeric variables 'SubjectID' and 'Generation'.
#' ExtraOutcomes79$SubjectTag <- CreateSubjectTag(
#'    subjectID=ExtraOutcomes79$SubjectID, 
#'    generation=ExtraOutcomes79$Generation
#' )
#' 
CreateSubjectTag <- function( subjectID, generation ) {
  if( base::length(subjectID) != base::length(generation) ) 
    base::stop("The length of the 'subjectID' vector did not match the length of the 'generation' vector.")
  
  tag <- base::rep(NA, base::length(subjectID))
   for( i in base::seq(base::length(subjectID)) ) {
    if( base::is.na(subjectID[i]) || base::is.na(generation[i]) )
      tag[i] <- NA
    else if( generation[i] == 1 ) 
      tag[i] <- subjectID[i] * 100L
    else if( generation[i] == 2 )
      tag[i] <- subjectID[i]
    else
      base::stop(base::paste("The generation value of '", generation[i], "' at element '", i, "' is not valid.  It must be either 1 or 2."))
   }
   return( tag )
}

# IncludeSubjectTag <- function( ds ) {
#   if( !("SubjectID" %in% colnames(ds)) ) stop("The data frame must contain a column named 'SubjectID' (case-sensitive).")
#   if( !("Generation" %in% colnames(ds)) ) stop("The data frame must contain a column named 'Generation' (case-sensitive).")  
#   ds$SubjectTag <<- CreateSubjectTag(subjectID=ds[, "SubjectID"], generation=ds[, "Generation"])
# }