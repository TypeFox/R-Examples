#' Get names of selections parameters from an summary.ezsim object.
#' @name getSelectionName.summary.ezsim
#' @aliases getSelectionName.summary.ezsim
#' @title Get Names of selection Parameters.
#' @usage \method{getSelectionName}{summary.ezsim}(x,simple=FALSE,parameters_priority,...)
#' @param x an summary.ezsim object
#' @param simple If true, return only the name of selection parameters. If False, split the selection into two groups, one with fixed value, one with varied value. Also, subtitle is returned.
#' @param parameters_priority Priority in sorting parameters.
#' @param \dots unused
#' @return Names of selection parameters. 
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export 
#' @note For internal use of ezsim.
#' @seealso \code{\link{getSelectionName.ezsim}}
#' keywords internal
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

getSelectionName.summary.ezsim<-
function(x,simple=FALSE,parameters_priority,...){
    selection_name <- attr(x,'selection_parameters')
    display_name <- attr(x,'display_name')
    
    selection_unique_length<-sapply(x[selection_name], function(x) length(unique(x)))
    
    selection_length_one<-names(selection_unique_length)[selection_unique_length==1]
    selection_length_greater_one<-selection_unique_length[selection_unique_length>1]
    
    selection_length_greater_one<-names(selection_length_greater_one)[order(selection_length_greater_one,decreasing=TRUE) ]
    
    subtitle<-''
    ## For selection_name_length_one, we update the subtitle.
    if (length(selection_length_one)>0 ){
        selection_value_length_one<-unique(x[selection_length_one])
        
        subtitle<-paste(paste(Jmisc::recode(selection_length_one,from=names(display_name),to=display_name),selection_value_length_one,sep='=='),collapse=',')
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
