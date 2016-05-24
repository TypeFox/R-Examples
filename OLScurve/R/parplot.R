#' Plot distribution of parameters
#' 
#' A plotting function for displaying the distribution of the OLS parameter
#' estimates. 
#' 
#' 
#' @aliases parplot
#' @param object an object of class \code{OLScurve}
#' @param type type of plot to display; can be \code{'hist'}, \code{'boxplot'}, or \code{'splom'}
#'    for a histogram, boxplot, or scatter plot matrix
#' @param group a \code{factor} grouping variable used to partition the results
#' @param breaks number of breaks to be used in plotting the histogram
#' @param prompt a logical variable indicating whether \code{devAskNewPage(ask=TRUE)} should be called
#' @param ... additional arguments to be passed
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords OLS, growth
#' @export parplot
#' @examples 
#' 
#' \dontrun{
#' data <- t(t(matrix(rnorm(1000),200)) + 1:5)  
#' group <- rep(c('Male', 'Female'), each=nrow(data)/2)
#' mod <- OLScurve(~ time, data = data)	
#' parplot(mod)
#' parplot(mod, type = 'boxplot')
#' parplot(mod, type = 'splom')
#' 
#' parplot(mod, group=group)
#' parplot(mod, type='boxplot', group=group)
#' parplot(mod, type='splom', group=group)
#' }
parplot <- function(object, ...){
	UseMethod('parplot')
}

#' @S3method parplot OLScurve
#' @rdname parplot 
#' @method parplot OLScurve 
parplot.OLScurve <- function(object, type = 'hist', group = NULL, 
	breaks = NULL, prompt = TRUE, ...)
{       
	pars <- as.data.frame(object$pars)    
    longpars <- data.frame(pars = as.numeric(object$pars),
                           coef = rep(colnames(object$pars), each=nrow(pars)))
    if(!is.null(group)){
        pars$group <- group
        longpars$group <- rep(group, ncol(object$pars))
    }
	if(type == 'splom'){			
        pars2 <- pars[, colnames(pars) != 'group']        
		if(is.null(group)) return(splom(~pars, data = pars, main = 'Growth Parameters'))
		else return(splom(~pars2|group, data = pars, main = 'Growth Parameters'))
	}	
	if(type == 'hist'){
	    if(is.null(group)) return(histogram(~pars|coef, data=longpars, breaks = breaks, 
                         main = 'Parameter Distributions'))			
	    else return(histogram(~pars|coef+group, data=longpars, breaks = breaks, 
                              main = 'Parameter Distributions'))
	}
	if(type == 'boxplot'){
	    if(is.null(group)) return(bwplot(~pars|coef, data=longpars, main = 'Parameter Distributions'))			
	    else return(bwplot(~pars|coef+group, data=longpars, main = 'Parameter Distributions'))    		
	}
}
