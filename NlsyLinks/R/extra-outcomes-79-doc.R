#' @name ExtraOutcomes79
#' @docType data

#' @title Extra outcome variables in the NLSY79
#' 
#' @description This dataset is provided primarily to facilitate documentation examples.
#'  
#' @format A data frame with 11,495 observations on the following 6 variables. There is one row per subject.  
#' \describe{ 
#' \item{SubjectTag}{The ID value assigned by NLS to the first subject.  For Gen1 Subjects, this is their "CaseID" (ie, R00001.00).  For Gen2 subjects, this is their "CID" (ie, C00001.00).} 
#' \item{SubjectID}{The ID value assigned by NLS to the first subject.  For Gen1 Subjects, this is their "CaseID" (ie, R00001.00).  For Gen2 subjects, this is their "CID" (ie, C00001.00).}
#' \item{Generation}{The generation of the subject.  Values are either 1 or 2, representing Gen1 and Gen2.  Note that this variable is not a  \code{factor} (in constrast with data frames like
#'    \code{\link{Links79Pair}}).  This dataset is supposed to mimick the dataset provided by the researcher, which typically will not have been converted to a \code{factor}.} 
#' \item{HeightZGenderAge}{The subject's height, standardized by gender and age (see Details).} 
#' \item{WeightZGenderAge}{The subject's weight, standardized by gender and age (see Details).} 
#' \item{AfqtRescaled2006Gaussified}{Armed Forces Qualification Test Score (Gen1 only; see Details).} 
#' \item{Afi}{Self-reported age of first intercourse (Gen1 only; see Details).} 
#' \item{Afm}{Self-reported age of first menstration (Gen1 only; see Details).}
#' \item{MathStandardized}{ Standardized PIAT Math scores (Gen2 only; see Details).} 
#' }
#' @author Will Beasley
#' @source Gen1 information comes from the Summer 2013 release of the
#' \href{http://www.bls.gov/nls/nlsy79.htm}{NLSY79 sample}.  Gen2 information comes from the Summer 2013 release of the
#' \href{http://www.bls.gov/nls/nlsy79ch.htm}{NLSY79 Children and Young Adults sample}.  Data were extracted with the NLS Investigator
#' (\url{https://www.nlsinfo.org/investigator/}).
#' 
#' @details 
#' The \code{SubjectTag} variable uniquely identify subjects.  For Gen2
#' subjects, the SubjectTag is identical to their CID (ie, C00001.00 -the
#' SubjectID assigned in the NLSY79-Children files).  However for Gen1
#' subjects, the SubjectTag is their CaseID (ie, R00001.00), with "00"
#' appended.  This manipulation is necessary to identify subjects uniquely in
#' inter-generational datasets.  A Gen1 subject with an ID of 43 has a
#' \code{SubjectTag} of 4300.  The SubjectTags of her four children remain
#' 4301, 4302, 4303, and 4304.
#' 
#' For Gen2, an NLSY79 variable of \code{MathStandardized} is C05801.00.
#' 
#' \code{Afi} and \code{Afm}, values were simplified
#' (to one value per subject) by Kelly Meredith in Sept 2010.
#' 
#' The variables for height and weight were manipulated in R files available in a 
#' \href{https://github.com/LiveOak/NlsyLinksDetermination/tree/master/ForDistribution/Outcomes}{repository} available to the public.
#' Find the appropriate subfolder, and view the HTML report for more details.
#' 
#' @keywords datasets
#' 
#' @examples 
#' library(NlsyLinks) #Load the package into the current R session.
#' gen2Outcomes <- subset(ExtraOutcomes79, Generation==2) #Create a dataset of only Gen2 subjects.
#'                   
#' #plot(ExtraOutcomes79) #Uncomment to see a large scatterplot matrix.
#' summary(ExtraOutcomes79)
#' 
#' oldPar <- par(mfrow=c(3,2))
#' hist(ExtraOutcomes79$Generation)
#' hist(ExtraOutcomes79$MathStandardized)
#' hist(ExtraOutcomes79$HeightZGenderAge)
#' hist(ExtraOutcomes79$WeightZGenderAge)
#' hist(ExtraOutcomes79$Afi)
#' hist(ExtraOutcomes79$Afm)
#' par(oldPar)
#' 
NULL
