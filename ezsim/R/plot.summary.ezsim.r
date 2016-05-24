#' Plot the summary statistics for several estimators in the same plot. Summary statistics abd estimators are separated by colour and linetype.
#' The longest scalars parameter will be the x-variable of the plot. The rest of the scalars parameters will be become the facets of the plot (see \pkg{ggplot2}).Banker parameters will not be shown in the graph.
#' @name plot.summary.ezsim
#' @aliases plot.summary.ezsim
#' @title Plot an summary.ezsim Object
#' @method plot summary.ezsim
#' @param x An summary.ezsim Object
#' @param parameters_priority Display priority of parameter. Any missed parameters will be sorted by length.
#' @param ylab Label of y-axis
#' @param title Title of the plot
#' @param return_print If TRUE, return a list of ggplot2 object. If FALSE(default), all of the plot will be printed out.
#' @param pdf_option A list of option pass to \code{\link{pdf}}. If it is not missing, the plot will export to a pdf file
#' @param \dots unused
#' @return Optional: a ggplot2 object
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export 
#' @seealso \code{\link{ezsim}},\code{\link{summary.ezsim}}
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
#' ## Plot the summary ezsim
#' plot(summary(ezsim_basic,c("q25","q75")))
#' plot(summary(ezsim_basic,c("q25","q75"),subset=list(estimator='mean_hat')))
#' plot(summary(ezsim_basic,c("median"),subset=list(estimator='sd_mean_hat')))
#' }
plot.summary.ezsim<-
function(x,parameters_priority,ylab='Summary Statistics',title,pdf_option,return_print,...){

    display_name<-attr(x,'display_name')
         
    temp<-getSelectionName(x,parameters_priority=parameters_priority)
    ########### title
    title<-
    if (missing(title)){
        temp$subtitle
    }
    else {  
        if (temp$subtitle!='')
            paste(title,temp$subtitle,sep='~~')
        else 
            title
    }   

    x_var=head(temp$selection_length_greater_one,1)
    other=tail(temp$selection_length_greater_one,-1)    
            
    summ<-x

    summ<-melt(summ,id.vars=c('estimator',getSelectionName(x,TRUE)),variable_name='stat')
    
		my_facet<-
			if (length(other>0)) {
				facet_grid(createFormula(other), labeller = Jmisc::label_both_parsed_recode(display_name))
			}else{
				NULL
			}
		
    out<-
        if (length(unique(summ$stat))==1 & length(unique(summ$estimator))==1){
            ggplot(data=summ, aes_string(x=x_var,y='value'),)+geom_line()
        } else if (length(unique(summ$stat))==1 & length(unique(summ$estimator))>1){
            ggplot(data=summ, aes_string(color='estimator', x=x_var,y='value'))+scale_colour_discrete(name='Estimators',breaks=unique(summ$estimator),labels=parse(text=paste(unique(summ$stat),'of', unique(summ$estimator),sep='~~')))
            
            #+scale_linetype(name='Summary Statistics',breaks=unique(summ$stat),labels=parse(text=unique(summ$stat)))
        } else if (length(unique(summ$stat))>1 & length(unique(summ$estimator))==1){
            ggplot(data=summ, aes_string(color='stat', x=x_var,y='value'))+scale_colour_discrete(name='Summary Statistics',breaks=unique(summ$stat),labels=parse(text=paste(unique(summ$stat),'of', unique(summ$estimator),sep='~~')))
            
            # +scale_linetype(name='Estimators',breaks=unique(summ$estimator),labels=parse(text=unique(summ$estimator)))
        } else {
            ggplot(data=summ, aes_string(linetype='stat',color='estimator', x=x_var,y='value'))+scale_colour_discrete(name='Estimators',breaks=unique(summ$estimator),labels=parse(text=unique(summ$estimator)))+scale_linetype(name='Summary Statistics',breaks=unique(summ$stat),labels=parse(text=unique(summ$stat)))
        }
        
    out <- out+
            geom_line()+
            geom_point() + 
            my_facet+
            ylab(ylab)+
            xlab(parse(text=recode(x_var,from=names(display_name),to=display_name)))+
            labs(title=parse(text=title))+
            theme(legend.position='bottom', legend.direction='horizontal')
        
    if (!missing(pdf_option)){
        do.call(pdf,pdf_option)
        if (class(out)=='ggplot')
            print(out)
        else
            lapply(out,print )
        dev.off()
    } 
    if (return_print){
        return(out)
    } else{
        if (class(out)=='ggplot')
            print(out)
        else
            temp<-lapply(out,function(x) {dev.new(); print(x)} )
    }
}
