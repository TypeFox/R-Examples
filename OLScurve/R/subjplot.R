#' Plot individually estimated parameters
#' 
#' A plotting function for displaying the individuals trajectories and their 
#' modelled functional form. Useful for detecting aberrant individual trajectories.
#' 
#' 
#' @aliases subjplot
#' @param object an object of class \code{OLScurve}
#' @param layout a variable to be passed to \code{xyplot} to adjust the graphical layout
#' @param prompt a logical variable indicating whether \code{devAskNewPage(ask=TRUE)} should be called
#' @param ... additional arguments to be passed
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords OLS, growth
#' @export subjplot
#' @examples 
#' 
#' \dontrun{
#' data <- t(t(matrix(rnorm(1000),200)) + 1:5)  
#' mod <- OLScurve(~ time, data = data)	
#' subjplot(mod)
#' 
#' ##quadratic
#' data <- t(t(matrix(rnorm(1000),200)) + (0:4)^2)
#' mod2 <- OLScurve(~ time + I(time^2), data = data)
#' subjplot(mod2)
#' 
#' 
#' }
subjplot <- function(object, ...){
	UseMethod('subjplot')
}

#' @S3method subjplot OLScurve
#' @rdname subjplot 
#' @method subjplot OLScurve 
subjplot.OLScurve <- function(object, layout = NULL, prompt = TRUE, ...)
{    
	data <- object$data
	N <- nrow(data)
	fn <- fn1 <- object$formula
	id.o <- id <- as.numeric(rownames(data))
	data <- data.frame(data)	
    colours <- rep(c('red','blue','black','darkviolet','green','black','goldenrod'), length.out=N)
	ypred <- object$pred
    lower <- object$lower
	upper <- object$upper
	yprednames <- colnames(ypred) <- paste('pred', 1:ncol(ypred),sep='')
    lowernames <- colnames(lower) <- paste('lower', 1:ncol(lower),sep='')
	uppernames <- colnames(upper) <- paste('upper', 1:ncol(upper),sep='')
	
	if(is.null(layout)) layout <- c(ceiling(log(N)),ceiling(log(N)))
	
	plotOLScurve <- function(data, fn, group = NULL, layout = NULL, pred)     {
        
		if(prompt) devAskNewPage(ask=TRUE)		
		data <- data.frame(data) 
		if(is.null(data$id)) data$id.o<-1:nrow(data) 
			else data$id.o <- data$id 
		if(!is.null(group)) 
			data$group <- group
		if(is.null(data$group)) { 
			data <- data[order(data$id.o),] 
			plot.fn <- y ~ time
		} else { 
			data <- data[order(data$group,data$id.o),]
			data$id.o <- paste(data$group,data$id.o) 
			plot.fn <- y ~ time | group
		}
		ys <- colnames(data)[!(colnames(data)=="id")&
			!(colnames(data)=="group")&
			!(colnames(data)=="id.o")]        		
		data2 <- data.frame(data,ypred,lower,upper,colours=colours)        
		datalg <- reshape(data2, idvar="id",
		                  varying = list(ys,yprednames,lowernames,uppernames),
		                  v.names = c("y","ypred","lower","upper"), 
		                  times = c(1:length(ys)),
		                  direction="long")
        
		ch <- as.character(fn)
		ch[2] <- gsub("x","y",ch[2],fixed=TRUE)
		ch[3] <- gsub("time","x",ch[3],fixed=TRUE)
		fn1 <- paste(ch[2],ch[1],ch[3])

		###### PLOT INDIVIDUAL PARTICIPANTS ######		
		mypanel = function(x, y, subscripts, lower, upper, colours, pred, ...){
		    upper <- upper[subscripts]
		    lower <- lower[subscripts]		    
		    panel.polygon(c(x,rev(x)),c(upper,rev(lower)), col=gray(.9), border=NA, ...)
            
			panel.xyplot(x, y, col = colours[subscripts], ...)            			
			panel.lines(x, pred[subscripts])						
		}
        subjectPlots <- xyplot( y ~ time | factor(id), data=datalg, 
                layout = layout,
    			xlab = "time",
    			ylab = "",
    			main = 'Subject plots',
                upper = datalg$upper,
    		    lower = datalg$lower,                  
                colours = datalg$colours, 
                pred=datalg$ypred,
    			panel = mypanel)                		  
		print(subjectPlots)        
	}	
	plotOLScurve(data, fn, group=NULL, layout, object$pred)	
	devAskNewPage(ask=FALSE)
}
