#' Return a subset of the simulation result of an ezsim object.
#' @name subset.ezsim
#' @aliases subset.ezsim
#' @title Return of the Simulation 
#' @method subset ezsim
#' @param x An ezsim Object
#' @param subset A list contains the subset of estimators and parameters. To select a subset of estimators: \code{list(estimator=c('name of estimator1','name of estimator2'))}. To select a subset of parameters: \code{list(mean=1:3,sd=4:5)}. Or both.
#' @param \dots unused
#' @return sim of ezsim
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @seealso \code{\link{ezsim}}
#' @note For internal use of ezsim.
#' @export 
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
#' subset(ezsim_basic,subset=list(estimator='mean_hat',mu=0,n=c(20,40)))
#' }
  
subset.ezsim<-
function(x,subset,...){
    if (is.null(x$simulation_table))
        stop('Please run the simulation first.')

    if (!missing(subset)){
        x$simulation_table<-
        tryCatch({
            if (length(subset)==0)
                stop('length of subset is zero \n')
            if (any(names(subset) %in% c('value_of_TV','value_of_estimator')))
                stop('subset cant contain  \'value_of_TV\' and \'value_of_estimator\' \n')
            if (!all(names(subset) %in% names(x$simulation_table) ))
                stop('Unknown name: ', names(subset)[!names(subset) %in% names(x$simulation_table)] ,'\n' )
                
            index<-apply(mapply(function(name,value) as.matrix(x$simulation_table[[name]] %in% value)  , name=names(subset) , value=subset),1, all)
            
            subset(x$simulation_table,subset=index)
        },
        error = function(e) {
            print(e)
            cat('Error in subseting. Nothing is dropped out. \n')
            x$simulation_table
        })
    }
    x
}
