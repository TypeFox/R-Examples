#'@title
#'Query Optimizely  to list experiments in projects
#'@description 
#'\code{GetExperimentList} extracts the experiments and summary information in a list of projects
#'@export 
#'
#'
#'@param project.list List of project identifiers
#'
#'
#'@examples
#'\dontrun{
#' # list experiments in a project 
#' # 
#' experiment.list<-GetExperimentList( c('1234567','7654321') )
#'}
#'
#'@return  data frame of all experiments in list of projects provided
#'
#'
#'@seealso \code{\link{GetExperimentInfo}}
#'
GetExperimentList <- function(project.list) {
    if (length(project.list) != 0) {
        experiment.df <- do.call(rbind, lapply(project.list, GetExperiment))
    }
    return(experiment.df)
} 
