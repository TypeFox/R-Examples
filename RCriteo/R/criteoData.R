#' @title Pull Criteo Campaign Report Data
#' 
#' @description This function manages the complete data download process. The requested data will be returned as a data frame.
#' Wrapper for \code{link{jobStatus}}, \code{link{getCriteoDownloadURL}} and \code{link{getCriteoData}}.
#' 
#' @param authToken Object created by \code{link{doCriteoAuth}}
#' @param appToken App Token as character string.
#' @param jobID Object created by \code{link{scedCriteoReport}}
#' 
#' @return Data
#' 
#' @export
criteoData <- function(authToken, appToken, jobID){
  #Simplyfies data download
  #wraps jobStatus(), getCriteoDownloadURL and getCriteoData
  #
  #arguments: authToken
  #           appToken
  #           jobID
  #
  #returns: dataframe
  jobStatus <- getCriteoJobStatus(authToken=authToken, appToken=appToken, jobID=jobID)
  while(jobStatus != "Completed"){ 
    jobStatus <- getCriteoJobStatus(authToken=authToken, appToken=appToken, jobID=jobID)
  }
  URL <- getCriteoDownloadURL(authToken = authToken, appToken = appToken, jobID = jobID)
  data <- getCriteoData(URL = URL, jobID = jobID)
  data
  
}