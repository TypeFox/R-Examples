#' @title Test information function
#' @export tff
#' @description plotting function for plotting the test information function (TIF).
#' @details no details in the moment.
#' @param pair_obj an object of class \code{"pair"} as a result from function \code{\link{pair}}.
#' @param items optional a vector (character or numeric) identifying the items (according their order in the data) to use for plotting the test information function.  
#' @param x The value(s) of the latent variable, at which the TIF will be evaluated. \code{x} should be either a numeric vector of theta values or a single numeric value. If \code{x} is given as a single numeric value plotting is supressed. If not given (default), 99 values spaced evenly between -4 and +4 will be used, handy for plotting.
#' @param main see parameters for \code{\link{plot}}
#' @param plot a logical (default \code{plot = TRUE}), defining wether to supress plotting an just return a matrix of the values of the Item information function.
#' @param cat a logical (default \code{cat = FALSE}), defining wether to plot as an overlay to the Test information function the item category information functions based on item categories. If \code{cat = TRUE} and \code{plot = FALSE} the values of the item category information functions are returned. 
#' @param lwd see parameters for \code{\link{plot}} 
#' @param col see parameters for \code{\link{plot}} 
#' @param ... arguments passed to plot
#' @return a plot, a "data.frame" or a single numeric with values of the Test information function.
#' @examples ########
#' data(sim200x3)
#' result <- pair(sim200x3)
#' tff(pair_obj = result) # TIF plot 
#' tff(pair_obj = result, cat=TRUE) # TIF plot 
#' tff(pair_obj = result, items=c("V1","V3"), cat=TRUE) # TIF plot 
#' tff(pair_obj = result, x=0) # TIF at theta=0 
#' tff(pair_obj = result, x=seq(0,4,.1)) # TIF for a given range of Thetas
#' ##### examples with other data ...
#' data(bfiN)
#' result <- pair(bfiN)
#' tff(pair_obj = result)
#' tff(pair_obj = result, cat=TRUE) # TIF plot 
####################################################

####################################################

tff <- function(pair_obj, items=NULL, x=NULL, main = "Test Information Function" , plot=TRUE, cat=FALSE, lwd = 2, col=1, ...  ){

# start 'ploting' function ------------   
if (length(items)==0){
  items <- 1:dim(pair_obj$threshold)[1]
}

itemname <- rownames(pair_obj$threshold[items,])

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

pvList <- list()
for (i in 1:length(items)){
  thres <- na.omit(pair_obj$threshold[items[i],])
  pvList[[i]] <- pvx.matrix(theta_v=theta_v, thres=thres)    
}
names(pvList) <- itemname

PV <- (data.frame(lapply(pvList,t)))

IFerg <- rowSums((PV*(1-PV)))
IFerg_cat <- (PV*(1-PV))

if (plot==TRUE){ 
  # leerplot
  plot(y= c(0, max(IFerg)) , x= c(min(theta_v),max(theta_v)) , ylim=c(0, max(IFerg)) ,bty="n",type="n", ylab=("Information"), xlab="Logits", main=main)
  ##### plotting -------------------------
  matlines(x=theta_v ,y=IFerg, lwd=lwd, col=col, ...)# ,...
  if (cat==TRUE){  
    matlines(x=theta_v ,y=IFerg_cat, lwd=lwd, col=col,...)# ,...
  }
  cat("Plotted Test Information Function", "\n", "with item(s): ", paste(itemname,collapse=", ") ,sep="")
}

# merke : 
#matplot(rowSums(t(pvx.matrix(theta_v=seq(-2,2,.1),thres=c(-.4,.1,.6))*(1-pvx.matrix(theta_v=seq(-2,2,.1),thres=c(-.4,.1,.6))))))

##### silent return -------------------------
if (plot==FALSE){(return(IFerg_cat))}
}


