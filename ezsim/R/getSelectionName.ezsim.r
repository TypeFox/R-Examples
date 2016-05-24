#' Get names of selection parameters from an ezsim object.
#' @name getSelectionName.ezsim
#' @aliases getSelectionName.ezsim
#' @title Get Names of selection Parameters.
#' @usage \method{getSelectionName}{ezsim}(x,simple=FALSE,parameters_priority,...)
#' @param x an ezsim object
#' @param simple If true, return only the name of selection parameters. If False, split the selection into two groups, one with fixed value, one with varied value. Also, subtitle is returned.
#' @param parameters_priority Priority in sorting parameters.
#' @param \dots unused
#' @return \item{selection_length_greater_one}{ Name of selection parameters with more than one elements} \item{selection_length_one}{Name of selection parameters with only one element} \item{subtitle}{subtitle for fixed selection parameters}
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @keywords internal
#' @note For internal use of ezsim.
#' @seealso \code{\link{getSelectionName.summary.ezsim}}
#' @examples       
#' \dontrun{
#' ezsim_basic<-ezsim(
#'     m             = 100,
#'     run           = TRUE,
#'     display_name  = c(mean_hat="hat(mu)",sd_mean_hat="hat(sigma[hat(mu)])"),
#'     parameter_def = createParDef(list(n=seq(20,80,20),mu=c(0,2),sigma=c(1,3,5))),
#'     dgp           = function() rnorm(n,mu,sigma),
#'     estimator     = function(x) c(mean_hat = mean(x), 
#'                                  sd_mean_hat=sd(x)/sqrt(length(x)-1)),
#'     true_value    = function() c(mu, sigma / sqrt(n-1))
#' )
#' 
#' getSelectionName(ezsim_basic)
#' getSelectionName(summary(ezsim_basic))
#' }

getSelectionName.ezsim <-
function(x,simple=FALSE,parameters_priority,...){
    selection_value<-unique(x$simulation_table[!names(x$simulation_table) %in% c('estimator','value_of_estimator','value_of_TV')])
    
    selection_unique_length<-sapply(selection_value, function(x) length(unique(x)))
    
    selection_length_one<-names(selection_value)[selection_unique_length==1]
    
    # sort rest of them
    selection_length_greater_one<-selection_unique_length[selection_unique_length>1]
    selection_length_greater_one<-names(selection_length_greater_one)[order(selection_length_greater_one,decreasing=TRUE) ]

    subtitle<-''
    ## For selections_name_length_one, we update the subtitle.
    if (length(selection_length_one)>0 ){
        selections_value_length_one<-unique(x$simulation_table[selection_length_one])
        
        subtitle<-paste(paste(Jmisc::recode(selection_length_one,from=names(x$display_name),to=x$display_name),selections_value_length_one,sep='=='),collapse=',')
        subtitle<-paste('list(',subtitle,')',sep='')
        
       subtitle<-paste('group(\'(\',',subtitle,',\')\' )')
    }
    
    ## Update selection_name_length_greater_one with parameters_priority
    if (!missing(parameters_priority)){
        parameters_priority<-
        tryCatch( match.arg(parameters_priority,selection_length_greater_one,TRUE) ,
            error = function(e) selection_length_greater_one )
        
        selection_length_greater_one<-c(parameters_priority,setdiff(selection_length_greater_one,parameters_priority))
    }
    
    if (!simple)
        return(list(selection_length_greater_one=selection_length_greater_one,selection_length_one=selection_length_one,subtitle=subtitle))
    else
        return(c(selection_length_greater_one,selection_length_one))
}
