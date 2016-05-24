#' A quick summary to the simulation. Summary statistics included mean, true value (tv), bias, bias percentage (mean/tv-1), sd, rmse (root mean square error), min, q25 (first quarter), median, q75 (third quarter), max, p value of jb-test. See \code{\link{ezsim}} and \code{\link{plot.summary.ezsim}} for details and examples.
#' @name summary.ezsim
#' @aliases summary.ezsim
#' @title Summarize an ezsim Object
#' @usage  
#' \method{summary}{ezsim}(object,stat=c('mean','tv','bias',
#' 'biaspercentage','sd','rmse','min','q25','median',
#' 'q75','max','jb_test'),simple=TRUE,subset,...)
#' @method summary ezsim
#' @param object An ezsim object
#' @param stat Some preset summary statistics. Included,  \code{c('mean','tv','bias','biaspercentage','sd','rmse','min','q25','median','q75','max','jb_test')}
#' @param simple If True, shows only mean, true value, bias, sd and rmse of the estimator. If False, shows all statistics in stat.
#' @param subset subset of estimators or parameters. See \code{\link{subset.ezsim}} for details.
#' @param \dots Furhter summary statistics. Given in the form stat_name=stat. For example, Mean=mean
#' @return A summary.ezsim object
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export 
#' @seealso \code{\link{ezsim}}, \code{\link{plot.summary.ezsim}}, \code{\link{getSelectionName.summary.ezsim}}
#' @keywords post-simulation
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
#' 
#' 
#' ## Summary of an ezsim object
#' summary(ezsim_basic)
#' 
#' ## Summary of a subset of ezsim object
#' summary(ezsim_basic,subset=list(estimator='mean_hat',n=c(20,40),sigma=c(1,3)))
#' 
#' ## More Summary Statistics
#' summary(ezsim_basic,simple=FALSE,subset=list(estimator='mean_hat',n=c(20,40),sigma=c(1,3)))
#' 
#' ## Customize the Summary Statistics
#' summary(ezsim_basic,stat=c("q25","median","q75"),Q025=quantile(value_of_estimator,0.025),
#'   Q975=quantile(value_of_estimator,0.975),subset=list(estimator='mean_hat',n=c(20,40),sigma=c(1,3)))
#' }

summary.ezsim <-
function(object,stat=c('mean','tv','bias','biaspercentage',
    'sd','rmse','min','q25','median','q75','max','jb_test'),simple=TRUE,subset,...){
    
    if (is.null(object$simulation_table))
        stop('Please run the simulation first.')
    if (length(stat)==0)
        stop('stat cant be empty')

    simulation_table<-subset.ezsim(object,subset=subset)$simulation_table

    ## preset summary statistics
    
    variable_name <- c('estimator',getSelectionName(object,TRUE)  )
		
    out<-NULL
    i<-NULL
    
    stat<-match.arg(stat, several.ok = TRUE)
    
    if (simple & length(stat)==12)
        stat=c('mean','tv','bias','sd','rmse','biaspercentage')

    if (length(stat)>0){
        stat_list<-
        foreach (i=stat,.combine=(function(...) paste(...,sep=',')) ) %do% {
            switch(i,
                mean='Mean=mean(value_of_estimator)',
                tv='TV=mean(value_of_TV)',
                bias='Bias=mean(value_of_estimator)-mean(value_of_TV)',
                biaspercentage='BiasPercentage=mean(value_of_estimator)/mean(value_of_TV)-1',
                sd='SD=sd(value_of_estimator)',
                rmse='rmse=sqrt(sd(value_of_estimator)^2+(mean(value_of_estimator)-mean(value_of_TV))^2)',
                min='Min=min(value_of_estimator)',
                q25='Q25=quantile(value_of_estimator,0.25)',
                median='Median=median(value_of_estimator)',
                q75='Q75=quantile(value_of_estimator,0.75)',
                max='Max=max(value_of_estimator)',
                jb_test='JB_test=Jmisc::JBTest(value_of_estimator)'
            )
        }

        out<-eval(parse(text=
            paste("ddply(simulation_table,variable_name,summarize,",stat_list,')',sep='')
        ))
    }
    
    ## Custom summary statistics
    arg_list<-names(lapply(match.call()[-1],deparse))

    if (! all(arg_list %in% c('object','stat','subset','simple')  )){
        temp<-ddply(simulation_table,variable_name,summarize,...)
        if (!is.null(out))
            out<-merge(out,temp,by=variable_name, suffixes = c("_1","_2"))
        else
            out<-temp
    }

    ## rename the estimator
    out$estimator<-factor(recode(as.character(out$estimator),from=names(object$display_name),to=object$display_name))

    ## Setup the summary.ezsim object
    # attr(out,'other_parameters')<-object$other_parameters
    attr(out,'selection_parameters')<-getSelectionName(object,TRUE)
    attr(out,'display_name')<-object$display_name

    class(out)<-c('summary.ezsim',class(out))

    out
}    
