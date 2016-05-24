#' There are 3 different modes to plot an ezsim object. \code{'summary'}, \code{'density'} and \code{'powerfun'} plot the summary statistics,density function and power function of an ezsim object respectively.\cr\cr
#' \code{'summary'}: The y-variable of the plot are summary statistics of the estimator. Two confidence bounds will be shaded in the plot. 25\% and 75\% percentile will form a 50\% confidence bound. Similarly, 2.5\% and 97.5\% percentile will form a 95\% confidence bound.  Each plot have only one estimator. The scalars parameter has the longest length will be the x-variable of the plot. The rest of the selection parameters will be become the facets of the plot (see \pkg{ggplot2}). \cr\cr \code{density} : Density plot of the estimator. Each plot have only one estimator. selection parameter will appear as different colour and in different facets.\cr\cr \code{powerfun} : Plot the power function of test(s). Estimators have to be a test (value = 1 if rejecting the null hypothesis, value = 0 if fail to reject the null hypothesis) banker parameters will not be shown in the graph.
#' @name plot.ezsim
#' @aliases plot.ezsim
#' @title Plot an ezsim Object
#' @method plot ezsim
#' @param x An ezsim object
#' @param type Type of plot
#' @param subset subset of estimators or parameters. See \code{\link{subset.ezsim}} for details.
#' @param parameters_priority Display priority of parameter. If any parameter is missing here, they will be sorted by length.
#' @param return_print If TRUE, return a list of ggplot2 object. If FALSE(default), all of the plot will be printed out.
#' @param ylab Label of y-axis
#' @param title Title of the plot
#' @param pdf_option A list of option pass to \code{\link{pdf}}. If it is not missing, the plot will export to a pdf file
#' @param null_hypothesis Null hypothesis of the test. For \code{type=='powerfun'} only.
#' @param benchmark Benchmark distribution. For \code{type=='density'} only.
#' @param \dots unused
#' @return Optional: a list of ggplot2 object
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @seealso \code{\link{ezsim}},\code{\link{summary.ezsim}}, \code{\link{plot.summary.ezsim}},
#' @examples       
#' \dontrun{
#' ## example 1
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
#' ## Plot an ezsim object
#' plot(ezsim_basic)
#' ## Subet of the Plot
#' plot(ezsim_basic,subset=list(estimator="sd_mean_hat",mu=0))
#' plot(ezsim_basic,subset=list(estimator="mean_hat",sigma=3))
#' ## Parameters Priority of the Plot
#' plot(ezsim_basic,subset=list(estimator="sd_mean_hat",mu=0),parameters_priority=c("sigma","n"))
#' plot(ezsim_basic,subset=list(estimator="mean_hat",sigma=c(1,3)),parameters_priority="mu")
#' 
#' ## Density Plot
#' plot(ezsim_basic,'density')
#' plot(ezsim_basic,"density",subset=list(estimator="mean_hat",sigma=3),parameters_priority="n",
#'    benchmark=dnorm)
#' plot(ezsim_basic,"density",subset=list(estimator="mean_hat",mu=0),parameters_priority="n" ,
#'    benchmark=dnorm)
#'  
#' ## example 2
#' ezsim_ols<-ezsim(
#'     m             = 100,    
#'     run           = TRUE,
#'     display_name  = c(beta_hat='hat(beta)',es='sigma[e]^2',xs='sigma[x]^2',
#'                       sd_beta_hat='hat(sigma)[hat(beta)]'),
#'     parameter_def = createParDef(selection=list(xs=c(1,3),beta=c(0,2),n=seq(20,80,20),es=c(1,3))),
#'     dgp           = function(){
#'                         x<-rnorm(n,0,xs)
#'                         e<-rnorm(n,0,es)
#'                         y<-beta * x + e
#'                         data.frame(y,x)
#'                     },
#'     estimator     = function(d){
#'                         r<-summary(lm(y~x-1,data=d))
#'                         out<-r$coef[1,1:2]
#'                         names(out)<-c('beta_hat','sd_beta_hat')
#'                         out
#'                     },
#'     true_value    = function() c(beta, es/sqrt(n)/xs) 
#' )
#' plot(ezsim_ols)
#' plot(ezsim_ols,subset=list(beta=0))
#' 
#' plot(ezsim_ols,'density')
#' plot(ezsim_ols,'density',subset=list(es=1,xs=1))
#' 
#'  
#' ## example 3
#' ezsim_powerfun<-ezsim(
#'     run           = TRUE,   
#'     m             = 100,
#'     parameter_def = createParDef(selection=list(xs=1,n=50,es=c(1,5),b=seq(-1,1,0.1))),
#'     display_name  = c(b='beta',es='sigma[e]^2',xs='sigma[x]^2'),
#'     dgp           = function(){
#'                         x<-rnorm(n,0,xs)
#'                         e<-rnorm(n,0,es)
#'                         y<-b * x + e
#'                         data.frame(y,x)
#'                     },
#'     estimator     = function(d){
#'                         r<-summary(lm(y~x-1,data=d))
#'                         stat<-r$coef[,1]/r$coef[,2]
#' 
#'                         # test whether b > 0
#'                         # level of significance : 5%
#'                         out <- stat > c(qnorm(.95), qt(0.95,df=r$df[2]))
#'                         names(out)<-c("z-test","t-test")
#'                         out
#'                     }
#' )
#' plot(ezsim_powerfun,'powerfun') 
#' }
plot.ezsim <-
function(x,type=c('summary','density','powerfun' ),subset,parameters_priority,return_print=FALSE,ylab,title,pdf_option,null_hypothesis,benchmark,...){

    if (is.null(x$simulation_table))
        stop('Please run the simulation first.')

    type=match.arg(type)
    out<-NULL
    i<-j<-Q025<-Q975<-Q25<-Q75<-TV<-Mean<-Median<-value_of_estimator<-NULL
    
    # selection_name<-getSelectionName(x)

    ##Summmary 
    if (type=='summary'){
        if (missing(ylab))
            ylab<-'Summary Statistics'
        
        ## Compute summary statistics
        summ<-summary(x,c('mean','q25','q75','tv','median'), Q025=quantile(value_of_estimator,0.025),simple=FALSE, Q975=quantile(value_of_estimator,0.975),subset=subset )
        
        temp<-getSelectionName(summ,parameters_priority=parameters_priority)
                
        if (missing(title))
            title=''
        subtitle<-temp$subtitle
            
        x_var=head(temp$selection_length_greater_one,1)
        other=tail(temp$selection_length_greater_one,-1)    
        
        #########
				my_facet<-
        if (length(other>0)) {
					facet_grid(createFormula(other), labeller = Jmisc::label_both_parsed_recode(x$display_name))
				}else{
					NULL
				}
				
        out<-dlply(summ,'estimator', 
            function(mydata) {
                mytitle<-''
                if (title=='')
                    mytitle<-paste('Summary','of',as.character(mydata$estimator),sep='~~')
                else
                    mytitle<-title
                if (subtitle!='')
                    mytitle<- paste(mytitle,subtitle,sep='~~')

                ggplot(aes_string(x=x_var),data=mydata)+

                geom_ribbon(aes(ymin=Q025,ymax=Q975,alpha=.1))+
                geom_ribbon(aes(ymin=Q25,ymax=Q75,alpha=.3),)+
                scale_alpha_continuous(name='Confidence Bound',range=c(.1,.3),breaks=c(0.1,.3),label=c('95%','50%'),limit=c(.1,.3))+
                
                geom_line(aes(y=TV,color='True Value'))+
                geom_line(aes(y=Mean,color='Mean'))+
                geom_line(aes(y=Median,color='Median'))+

                geom_point(aes(y=TV,color='True Value'))+
                geom_point(aes(y=Mean,color='Mean'))+
                geom_point(aes(y=Median,color='Median'))+
                
                scale_colour_manual(name='Summary Statistics',values=c('red','black','blue'))+
                my_facet + ylab(ylab)+xlab(parse(text=Jmisc::recode(x_var,from=names(x$display_name),to=x$display_name)))+
                theme(legend.position='bottom', legend.direction='horizontal') +
                labs(title = parse(text=mytitle))     
            }
        )
    }
    
    ## powerfun
    if (type=='powerfun'){
        if (missing(ylab))
            ylab<-'Probability of Rejecting Null Hypothesis'
            summary_x<- summary(x,stat='mean',subset=subset,simple=FALSE)


        mytitle<-
        if (missing(title))
            'Power~~Function'
        else
            title          

        out<-plot(summary_x,parameters_priority=parameters_priority,return_print=TRUE,ylab=ylab ,title=mytitle)
        
       
        if (!missing(null_hypothesis)){
            out<-out+
            eval(substitute(geom_vline(xintercept=null),list(null=null_hypothesis)))
        }
    }
    
    ## density
    if (type=='density'){
        if (missing(ylab))
            ylab<-'Density'
        
        x$simulation_table<-subset.ezsim(x,subset=subset)$simulation_table
        
        selection_name <- getSelectionName(x,parameters_priority=parameters_priority)
 
        x_var=head(selection_name$selection_length_greater_one,1)
        other=tail(selection_name$selection_length_greater_one,-1)    

        # change name of estimator to display name
        x$simulation_table$estimator<-factor(Jmisc::recode(as.character(x$simulation_table$estimator),from=names(x$display_name),to=x$display_name))
        
				my_facet<-
					if (length(other>0)) {
						facet_grid(createFormula(other,right=FALSE), labeller = Jmisc::label_both_parsed_recode(x$display_name))
					}else{
						NULL
					}
        
        if (!missing(benchmark) ){
            benchmark<-
            if (length(benchmark)!= length(unique(x$simulation_table$estimator)))
                replicate(length(unique(x$simulation_table$estimator)),benchmark)
            else
                list(benchmark)
        }    
            
        out<-
        foreach ( i = unique(x$simulation_table$estimator), j=1:length(unique(x$simulation_table$estimator)) ) %do% {
            temp<-subset.data.frame(x$simulation_table,subset=x$simulation_table$estimator==i)
            temp[[x_var]]<-factor(temp[[x_var]])
                                
            mytitle<-as.character(i)
            
            if (missing(title))
                mytitle<- paste('Density~~of',mytitle,sep='~~')
            else
                mytitle<- paste(title,'Density~~of',mytitle,sep='~~')
                
            if (selection_name$subtitle!='')
                mytitle<-paste(mytitle,selection_name$subtitle,sep='~~')

            temp_out<-
            ggplot(data=temp)+
            geom_density(aes_string(x='value_of_estimator',color=x_var,fill=x_var),alpha=0.1) + 
            scale_color_discrete(guide="none") + 
            scale_fill_discrete(name=parse(text=x_var)) +
            my_facet + 
            ylab(ylab) + 
            xlab(parse(text=as.character(i))) +
            theme(legend.position='bottom', legend.direction='horizontal') +
            labs(title=parse(text=mytitle))
          
            
            
            if (!missing(benchmark)){
                temp_out<-temp_out+
                stat_function(fun = dnorm, colour="black",aes_string(x='value_of_estimator')) 
            }
            temp_out
        }
    }
        
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
