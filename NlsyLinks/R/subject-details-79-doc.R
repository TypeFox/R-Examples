#' @name SubjectDetails79
#' @docType data

#' @title Dataset containing further details of the Gen1 and Gen2 subjects.
#' 
#' @description These variables are useful to many types of analyses (not just behavior genetics), and are provided to save users time.
#' 
#' 
#' @name SubjectDetails79
#' @docType data
#' @format A data frame with 24,181 observations on the following 12
#' variables.  
#' \describe{ 
#'    \item{SubjectTag}{see the variable of the same name in \code{\link{Links79Pair}}} 
#'    \item{ExtendedID}{see the variable of the same name in \code{\link{Links79Pair}}}
#'    \item{Generation}{Indicates if the subject is in generation \code{1} or \code{2}.} 
#'    \item{Gender}{Indicates if the subject is \code{Male} or \code{Female}.} 
#'    \item{RaceCohort}{Indicates if the race cohort is \code{Hispanic}, \code{Black} or \code{Nbnh} (\emph{ie}, Non-black, non-hispanic).  This comes from the Gen1 variable \code{R02147.00} and Gen2 variable \code{C00053.00}.} 
#'    \item{SiblingCountInNls}{The number of the subject's siblings, including himself/herself (a singleton has a value of one).  This considers only the siblings in the NLSY.  For Gen1, this can exclude anyone outside the age range.  For Gen2, this excludes anyone who doesn't share the same mother.} 
#'    \item{BirthOrderInNls}{Indicates the subject's birth order among the NLSY siblings.}
#'    \item{SimilarAgeCount}{The number of children who were born within roughly 30 days of the subject's birthday, including the subject (for instance, even an only child will have a value of 1).  For Gen2 subjects, this should reflect how many children the Gen1 mother gave birth to at the same time (1: singleton; 2: twins, 3: triplets).  For Gen1 subjects, this is less certain, because the individual might have been living with a similarly-aged housemate, born to a different mother.}
#'    \item{HasMzPossibly}{Indicates if the subject \emph{might} be a member of an MZ twin/triplet. This will be true if there is a sibling with a DOB within a month, and they are the same gender.}
#'    \item{IsMz}{Indicates if the subject has been identified as a member of an MZ twin/triplet. }
#'    \item{KidCountBio}{The number of biological children known to the NLSY (but not necessarily interviewed by the NLSY.}
#'    \item{KidCountInNls}{The number of children who belong to the NLSY. This is nonnull for only Gen1 subjects.} 
#'    \item{Mob}{The subject's month of birth.  The exact day is not available to the public.  By default, we set their birthday to the 15th day of the month.}
#'    \item{LastSurveyYearCompleted}{The year of the most recently completed survey.}
#'    \item{AgeAtLastSurvey}{The subject's age at the most recently completed survey.}
#'    \item{IsDead}{##This variable is not available yet## Indicates if the subject was alive for the last attempted survey.}
#'    \item{DeathDate}{##This variable is not available yet## The subject's month of death.  The exact day is not available to the public. By default, we set their birthday to the 15th day of the month.} 
#' }
#' 
#' @author Will Beasley
#' @source Gen1 information comes from the Summer 2013 release of the \href{http://www.bls.gov/nls/nlsy79.htm}{NLSY79 sample}.  Gen2 information
#' comes from the Summer 2013 release of the
#' \href{http://www.bls.gov/nls/nlsy79ch.htm}{NLSY79 Children and Young Adults sample}.  Data were extracted with the NLS Investigator
#' (\url{https://www.nlsinfo.org/investigator/}).
#' @keywords datasets
#' @examples 
#' library(NlsyLinks) #Load the package into the current R session.
#' 
#' summary(SubjectDetails79)
#' 
#' oldPar <- par(mfrow=c(3,2), mar=c(2,2,1,.5), tcl=0, mgp=c(1,0,0))
#' hist(SubjectDetails79$SiblingCountInNls, main="",
#'      breaks=seq(from=0, to=max(SubjectDetails79$SiblingCountInNls, na.rm=TRUE), by=1)
#' )
#' hist(SubjectDetails79$BirthOrderInNls, main="",
#'      breaks=seq(from=0, to=max(SubjectDetails79$BirthOrderInNls, na.rm=TRUE), by=1)
#' )
#' hist(SubjectDetails79$SimilarAgeCount, main="",
#'      breaks=seq(from=0, to=max(SubjectDetails79$SimilarAgeCount, na.rm=TRUE), by=1)
#' )
#' hist(SubjectDetails79$KidCountBio, main="",
#'      breaks=seq(from=0, to=max(SubjectDetails79$KidCountBio, na.rm=TRUE), by=1)
#' )
#' hist(SubjectDetails79$KidCountInNls, main="",
#'      breaks=seq(from=0, to=max(SubjectDetails79$KidCountInNls, na.rm=TRUE), by=1)
#' )
#' #hist(SubjectDetails79$Mob, main="", 
#' #     breaks=seq.Date(
#' #       from=min(SubjectDetails79$Mob, na.rm=TRUE), 
#' #       to=max(SubjectDetails79$Mob, na.rm=TRUE), 
#' #       by="year")
#' #)
#' par(oldPar)
#' 
NULL