#' @name RGroupSummary
#' @export
#' 
#' @title Calculates summary statistics for each \emph{R}elatedness Group in the sample.
#' 
#' @description Before and after running ACE Models, it is important to examine the characteristics of the different groups.  When the ACE is estimated with an SEM using multiple groups, it is even even more important.  Groups may contain too few subjects to have a well-behaved covariance matrix.  
#' 
#' If a group's covariance matrix is not Positive Definite (or it's misbehaving in some other way), it's typically recommended to exclude that group from the SEM. 
#' @usage RGroupSummary(ds, oName_S1, oName_S2, rName, determinantThreshold=1e-5)
#' 
#' @param ds The \code{data.frame} containing the following variables:
#' @param oName_S1 The name of the outcome variable corresponding to the first subject in the pair.
#' @param oName_S2 The name of the outcome variable corresponding to the first subject in the pair.
#' @param rName The name of the variable specifying the pair's \code{Relatedness} coefficient.
#' @param determinantThreshold The minimum value the covariance matrix's determinant (for the group) should exceed to be considered Positive Definite.
#' 
#' @details This function isn't specific to an ACE model and groups defined by \code{R}.  It could be applied to any multiple-group SEM with two manifest/outcome variables.  In the future, we may generalize it beyond two manifest variables.
#' 
#' To get summary stats for the entire sample, create a dummy indicator variable that assigns everyone to the same group.  See the second example below.
#' 
#' The default determinantThreshold value is nonzero, in order to forgive slight numerical inaccuracies caused by fixed-precision arithmetic.
#' 
#' @return A \code{data.frame} with one row per group.  The \code{data.frame} contains the following variables:
#' \item{ R }{ The group's \code{R} value.  Note the name of this variable can be changed by the user, by specifying a non-default value to the \code{rName} argument.}
#' \item{ Included }{ Indicates if the group should be included in a multiple-group SEM.}
#' \item{ PairCount }{ The number of pairs in the group with \emph{complete} data for \code{R} and the two outcome/manifest variables.}
#' \item{ O1Mean }{ The mean (of the outcome variable) among the group's first members, excluding the missing values. }
#' \item{ O2Mean }{ The mean (of the outcome variable) among the group's second members, excluding the missing values. }
#' \item{ O1Variance }{ The variance (of the outcome variable) among the group's first members. }
#' \item{ O2Variance }{ The variance (of the outcome variable) among the group's second members. }
#' \item{ O1O2Covariance }{ The covariance (of the outcome variable) across the group's first and second members.}
#' \item{ Correlation }{ The correlation (of the outcome variable) across the group's first and second members.}
#' \item{ Determinant }{ The determinant of the group's covariance matrix.}
#' \item{ PosDefinite }{ Indicates if the group's covariance matrix is positive definite.}
#' 
#' @references Please see \href{http://ibgwww.colorado.edu/workshop2006/cdrom/HTML/book2004a.pdf}{Neale & Maes} for more information about SEM with multiple groups.
#' 
#' @author Will Beasley and David Bard
#' 
#' @examples
#' library(NlsyLinks) #Load the package into the current R session.
#' dsLinks <- Links79PairExpanded  #Load the dataset from the NlsyLinks package.
#' dsLinks <- dsLinks[dsLinks$RelationshipPath=='Gen2Siblings', ]
#' oName_S1 <- "MathStandardized_S1" #Stands for Outcome1
#' oName_S2 <- "MathStandardized_S2" #Stands for Outcome2
#' dsGroupSummary <- RGroupSummary(dsLinks, oName_S1, oName_S2)
#' dsGroupSummary
#' 
#' #Should return: 
#' #      R Included PairCount   O1Mean   O2Mean O1Variance O2Variance O1O2Covariance Correlation
#' #1 0.250     TRUE      2718  94.6439  95.5990    169.650    207.842        41.0783    0.218761
#' #2 0.375     TRUE       139  92.6043  93.1655    172.531    187.081        40.4790    0.225311
#' #3 0.500     TRUE      5511  99.8940 100.1789    230.504    232.971       107.3707    0.463336
#' #4 0.750    FALSE         2 108.5000 106.0000    220.500     18.000        63.0000    1.000000
#' #5 1.000     TRUE        22  98.6364  95.5455    319.195    343.117       277.5887    0.838789
#' #  Determinant PosDefinite
#' #1     33573.0        TRUE
#' #2     30638.7        TRUE
#' #3     42172.2        TRUE
#' #4         0.0       FALSE
#' #5     32465.6        TRUE
#' 
#' #To get summary stats for the whole sample, create one large inclusive group.
#' dsLinks$Dummy <- 1
#' (dsSampleSummary <- RGroupSummary(dsLinks, oName_S1, oName_S2, rName="Dummy"))
#'                      
#' #Should return:
#' #  Dummy Included PairCount   O1Mean   O2Mean O1Variance O2Variance O1O2Covariance
#' #1     1     TRUE      8392 98.07162 98.56864    216.466   229.2988       90.90266
#' #  Correlation Determinant PosDefinite
#' #1   0.4080195     41372.1        TRUE
#' ###
#' ### ReadCsvNlsy79
#' ###
#' \dontrun{
#' filePathGen2 <- "~/Nlsy/Datasets/gen2-birth.csv"
#' ds <- ReadCsvNlsy79Gen2(filePath=filePathGen2)
#' }
#' 
#' @keywords ACE 

RGroupSummary <- function( ds, oName_S1, oName_S2, rName="R", determinantThreshold=1e-5) {
#     ds <- Links79PairExpanded #Start with the built-in data.frame in NlsyLinks
#     oName_S1 <- "MathStandardized_S1" #Stands for Manifest1
#     oName_S2 <- "MathStandardized_S2" #Stands for Manifest2
#   
#   ds <-dsFull
#   rName <- "RRR"
  
  #ds <- subset(ds, R==.75)
  rLevelsFirstPass <- base::sort(base::unique(ds[, rName])) #Enumerate the values of R existing in the current data.frame.
   dsGroupSummary <- base::data.frame(
    R=rLevelsFirstPass, Included=F, PairCount=NA, O1Mean=NA, O2Mean=NA,
    O1Variance=NA, O2Variance=NA, O1O2Covariance=NA, Correlation=NA, 
    Determinant=NA, PosDefinite=FALSE)
    
  index <- VerifyColumnExists(dataFrame=dsGroupSummary, columnName="R")
  base::colnames(dsGroupSummary)[index] <- rName
  
  #The primary goal of this loop is to identify the R groups whose covariance matrix isn't positive definite.
  #TODO: Consider rewriting with plyr.  It's small, so there won't be a performance benefit.  It will add another dependency to the package.
  for( rLevel in rLevelsFirstPass ) {
    #print(rLevel)
    dsGroupSlice <- ds[!is.na(ds[, rName]) & ds[, rName]==rLevel & !is.na(ds[, oName_S1]) & !is.na(ds[, oName_S2]), c(oName_S1, oName_S2)]
    
    if( base::nrow(dsGroupSlice) > 0 ) {
      o1Mean <- base::mean(dsGroupSlice[, oName_S1], na.rm=TRUE)
      o2Mean <- base::mean(dsGroupSlice[, oName_S2], na.rm=TRUE)
      groupCovarianceMatrix <- stats::cov(dsGroupSlice)#, use="complete.obs") 
      determinant <- base::det(groupCovarianceMatrix)
      isPositiveDefinite <- (determinant > determinantThreshold)
      correlation <- stats::cor(dsGroupSlice[, oName_S1], dsGroupSlice[, oName_S2])
    }
    else {
      o1Mean <- NA
      o2Mean <- NA
      groupCovarianceMatrix <- base::matrix(NA, ncol=2, nrow=2)
      determinant <- NA
      isPositiveDefinite <- F
      correlation <- NA
    }
    
    dsGroupSummary[dsGroupSummary[, rName]==rLevel, c("PairCount", "O1Mean", "O2Mean", "O1Variance", "O2Variance", "O1O2Covariance", "Correlation", "Determinant", "PosDefinite")] <- c(
      nrow(dsGroupSlice),
      o1Mean,
      o2Mean,
      groupCovarianceMatrix[1, 1],
      groupCovarianceMatrix[2, 2],
      groupCovarianceMatrix[1, 2],
      correlation,
      determinant,
      isPositiveDefinite
      )    
    #dsGroupSummary[dsGroupSummary[,rName]==rLevel, "Included"] <- isPositiveDefinite  
  }
  dsGroupSummary$PosDefinite <- base::as.logical(dsGroupSummary$PosDefinite) #I do not know how this variable was ever coerced from logical to numeric.
  dsGroupSummary[, "Included"] <- dsGroupSummary$PosDefinite #Maybe later there will be another criterion to include/exclude a group.
  
  return( dsGroupSummary )
}
