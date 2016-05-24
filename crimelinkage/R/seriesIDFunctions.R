
## seriesID
##==============================================================================
#' Crime series identification
#'
#' Performs crime series identification by finding the crime series that are
#'  most closely related (as measured by Bayes Factor) to an unsolved crime.
##  Inputs:
#'  @param crime crime incident; vector of crime variables
#'  @param solved incident data for the solved crimes. Must have a column named 
#'   \code{crimeID}.
#'  @param seriesData table of crimeIDs and crimeseries (results from 
#'  \code{\link{makeSeriesData}})
##  Case Linkage Objects:
#'  @param varlist a list of the variable names (columns of \code{solved} and 
#'   \code{crime}) 
#'   used to create evidence variables with \code{\link{compareCrimes}}.   
#'  @param estimateBF function to estimate the bayes factor from evidence variables
#'  @param linkage.method the type of linkage for comparing one crime to a set 
#'   of crimes
#'    \itemize{
#'      \item \dQuote{average} uses the average bayes factor
#'      \item \dQuote{single} uses the largest bayes factor (most similar)
#'      \item \dQuote{complete} uses the smallest bayes factor (least similar)
#'    }
#'  @param group.method the type of crime groups to form (see \code{\link{makeGroups}} 
#'    for details)
#'  @param \ldots other arguments passed to \code{\link{compareCrimes}}
##  Outputs:
#'  @return A list with two objects. \code{score} is a data.frame of the similarity 
#'    scores for each element in \code{solved}. \code{groups} is the data.frame
#'    \code{seriesData} with an additional column indicating the crime group 
#'    (using the method specified in \code{group.method}).
#'  @examples
#'  # See vignette: "Crime Series Identification and Clustering" for usage.    
#'  @references
#'  Porter, M. D. (2014). A Statistical Approach to Crime Linkage.
#'    \emph{arXiv preprint arXiv:1410.2285.}.
#'  \url{http://arxiv.org/abs/1410.2285}
#'  @export
##==============================================================================
seriesID <- function(crime,solved,seriesData,varlist,estimateBF,
                     linkage.method = c('average','single','complete'),
                     group.method = 3,...){

  #- make crime data
  vars = unique(c('crimeID',unlist(varlist)))
  if(is.null(crime$crimeID)) crime$crimeID = 'new.crime'
  crimedata = rbind(crime[,vars],solved[,vars])  # add new.crime to solved

  #- compare crime pairs, get evidence vars
  pairs = data.frame(i1=crime$crimeID,i2=unique(solved$crimeID),stringsAsFactors=FALSE) 
  X = compareCrimes(pairs,crimedata,varlist,...)  

  #- Estimate the Bayes Factor
  bf = estimateBF(X)  

  #- Find 'nearest' series to new.crime
  CG = makeGroups(seriesData,method=group.method)          
  SD = unique(data.frame(crimeID=seriesData$crimeID,CG))
  SD$BF = bf[match(SD$crimeID,pairs$i2)]
  SD = na.omit(SD)
  Y = linkage(SD$BF,SD$CG,method=linkage.method)
return(list(score=Y,groups=data.frame(seriesData,group=CG)))
}


## linkage
##==============================================================================
#' Hierarchical Based Linkage
#'
#' Groups the Bayes Factors by crime group and calculates the linkage score for 
#'   each group. 
##  Inputs:
#'  @param BF vector of Bayes Factors
#'  @param group crime group
#'  @param method the type of linkage for comparing a crime to a set of crimes
#'    \itemize{
#'      \item \dQuote{average} uses the average bayes factor
#'      \item \dQuote{single} uses the largest bayes factor (most similar)
#'      \item \dQuote{complete} uses the smallest bayes factor (least similar)
#'    }
##  Outputs:
#'  @details If methods is a vector of linkages to use, then the all linkages are 
#'    calcualted and ordered according to the first element.
#'  @return a data.frame of the Bayes Factor scores
#'    ordered (highest to lowest). 
#'  @examples
#'  # See vignette: "Crime Series Identification and Clustering" for usage.    
#'  @export
##==============================================================================
linkage <- function(BF,group,method=c('average','single','complete')){
  method = match.arg(method,several.ok=TRUE)
  link.fun <- function(x) c( single   = max(x,na.rm=TRUE),
                             complete = min(x,na.rm=TRUE),
                             average  = mean(x,na.rm=TRUE))
  Y = do.call(rbind,tapply(BF,group,link.fun))
  Y = data.frame(group=rownames(Y),Y[,method,drop=FALSE])
  Y = Y[order(Y[,method[1]],decreasing=TRUE),]    
  rownames(Y) <- NULL
return(Y)
}


# 
# 
# ## seriesID.batch
# ##==============================================================================
# #' Crime series identification (batch mode)
# #'
# #' Performs crime series identification by finding the crime series that are
# #'  most closely related to a \bold{set} of unsolved crimes.
# #'
# #'  Runs \code{\link{seriesID}} for all unsolved crimes
# ##  and return only the highest ranked crimes
# ##  Inputs:
# #'  @param unsolved set of unsolved crime incidents
# #'  @param solved crime incident data of solved crimes
# #'  @param seriesData table of crimeIDs and crimeseries (results from makeSeriesData)
# ##  Case Linkage Objects:
# #'  @param varlist the variable names necessary for getting evidence variables
# #'  @param binary (logical) match/no match or all combinations for categorical class
# #'  @param estimateBF function to estimate the bayes factor from evidence variables
# #'  @param linkage the type of linkage for compare a crime to a set of crimes
# #'    \itemize{
# #'      \item \dQuote{average} uses the average bayes factor
# #'      \item \dQuote{single} uses the largest bayes factor
# #'      \item \dQuote{complete} uses the smallest bayes factor
# #'    }
# #'  @param show.pb (logical) should a progress bar be displayed
# ##  Outputs:
# #'  @return a list where each element is the result of seriesID with
# #'    \code{details = FALSE}
# ##  @export
# ##==============================================================================
# seriesID.batch <- function(unsolved,solved,seriesData,varlist,estimateBF,
#                             linkage=c('average','single','complete'),
#                             binary=TRUE,show.pb=TRUE){
#   n = nrow(unsolved)
#   score = vector('list',n); names(score) = paste0('crime',1:n)
#   if(show.pb) pb = txtProgressBar(style=3,max=n)
#   for(i in 1:n){
#     score[[i]] = seriesID(unsolved[i,],solved,seriesData,varlist,estimateBF,
#                           linkage=linkage,details=FALSE)
#     if(show.pb) setTxtProgressBar(pb,i)
#   }
#   if(show.pb) close(pb)
# return(score)
# }
# 

