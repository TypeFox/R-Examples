#' @title Item information function
#' @export iff
#' @description plotting function for plotting the Item information function(IIF).
#' @details no details in the moment.
#' @param pair_obj an object of class \code{"pair"} as a result from function \code{\link{pair}}.
#' @param itemnumber an integer, defining the number of the item to plot the respective item information function for. This is set to an arbitrary default value of \code{itemnumber = 1} to avoid error messages when you forget to choose an item to plot the item information function for.
#' @param x The value(s) of the latent variable, at which the IIF will be evaluated. \code{x} should be either a numeric vector of theta values or a single numeric value. If \code{x} is given as a single numeric value plotting is supressed. If not given (default), 99 values spaced evenly between -4 and +4 will be used, handy for plotting.
#' @param plot a logical (default \code{plot = TRUE}), defining wether to supress plotting an just return a matrix of the values of the Item information function.
#' @param cat a logical (default \code{cat = FALSE}), defining wether to plot or return the values of the Item information function based on item categories.
#' @param lwd see parameters for \code{\link{plot}} 
#' @param col see parameters for \code{\link{plot}} 
#' @param ... arguments passed to plot
#' @return a plot, a matrix or a single numeric with values of the Item information function.
#' @examples ########
#' data(sim200x3)
#' result <- pair(sim200x3)
#' # IFF plot for Item No. 2 
#' iff(pair_obj = result, itemnumber = 2 ) 
#' # IFF plot for Categories of Item No. 2
#' iff(pair_obj = result, itemnumber = 2 ,cat=TRUE)
#' # IFF at theta=0 for Item No. 2 
#' iff(pair_obj = result, itemnumber = 2 ,x=0) 
#' # IFF at theta=0 for Categories of Item No. 2
#' iff(pair_obj = result, itemnumber = 2 ,x=0,cat=TRUE)
#' # IFF of Item No. 2 for a given range of thetas 
#' iff(pair_obj = result, itemnumber = 2 ,x=seq(0,4,.1)) 
#' # ... etc.
#' iff(pair_obj = result, itemnumber = 2 ,x=seq(0,4,.1),cat=TRUE) 
#' ##### examples with other data ...
#' data(bfiN)
#' result <- pair(bfiN)
#' iff(pair_obj = result, itemnumber = 3 )
#' iff(pair_obj = result, itemnumber = 3 ,cat=TRUE)
####################################################

####################################################

iff <- function(pair_obj, itemnumber=1, x=NULL, plot=TRUE, cat=FALSE, lwd = 2, col=1, ...  ){

# start 'ploting' function ------------    
itemname <- rownames(pair_obj$threshold)[itemnumber]

if (length(x)==0){ 
  ra <- 4
  theta_v <- seq(-ra,ra,length.out=(ra*2+1)*10)
}

if (length(x)==1){
  plot <- FALSE
  theta_v <- x
}

if (length(x)>1){
  theta_v <- x
}

thres <- na.omit(pair_obj$threshold[itemnumber,])
pv <- pvx.matrix(theta_v=theta_v, thres=thres)

if (cat==FALSE){IFerg <- rowSums(t(pv*(1-pv))) }
if (cat==TRUE){IFerg <- t(pv*(1-pv)) }

if (plot==TRUE){ 
  # leerplot
  plot(y= c(0, max(IFerg)) , x= c(min(theta_v),max(theta_v)) , ylim=c(0, max(IFerg)) ,bty="n",type="n", ylab=("Information"), xlab="Logits", main=paste("Item Information Function \n for item ", itemname,sep=""))
  ##### plotting -------------------------
  matlines(x=theta_v ,y=IFerg, lwd=lwd, col=col, ... ) 
}

# merke : 
#matplot(rowSums(t(pvx.matrix(theta_v=seq(-2,2,.1),thres=c(-.4,.1,.6))*(1-pvx.matrix(theta_v=seq(-2,2,.1),thres=c(-.4,.1,.6))))))

##### silent return -------------------------
if (plot==FALSE){(return(IFerg))}
}


