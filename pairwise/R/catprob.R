#' @title Category Probability Plots
#' @export catprob
#' @description plotting function for plotting category probability curves.
#' @details no details in the moment.
#' @param pair_obj an object of class \code{"pair"} as a result from function \code{\link{pair}}.
#' @param itemnumber an integer, defining the number of the item to plot the respective category probability for. This is set to an arbitrary default value of \code{itemnumber = 1} to avoid error messages when you forget to choose an item to plot the expected score curves for.
#' @param ra an integer, defining the (logit) range for x-axis
#' @param plot a logical (default \code{plot = TRUE}), defining wether to supress plotting an just return a matrix of category probabilities 
#' @param ... arguments passed to plot
#' @return a plot or a matrix with category probabilities.
#' @examples ########
#' data(sim200x3)
#' result <- pair(sim200x3)
#' catprob(pair_obj = result, itemnumber = 2 )
#' data(bfiN)
#' result <- pair(bfiN)
#' catprob(pair_obj = result, itemnumber = 3 )
####################################################


####################################################

catprob <- function(pair_obj, itemnumber=1, ra=4, plot=TRUE, ...  ){

# start 'ploting' function ------------    
itemname <- rownames(pair_obj$threshold)[itemnumber]

p <- round(t(sapply(seq(-ra,ra,length.out=(ra*2+1)*10), function(x){pvx(theta=x,thres=na.omit(pair_obj$threshold[itemnumber,]))} )),9)

rownames(p) <- seq(-ra,ra,length.out=(ra*2+1)*10)

dim(p) 
##### plotting -------------------------

if (plot==TRUE){ 
plot(y=p[,1], x= seq(-ra,ra,length.out=(ra*2+1)*10), ylim=c(0,1) ,bty="n",type="n",xaxt="n", ylab=("p"), xlab="logits", main=paste("category probability curves for item ", itemname,sep=""))
for (i in 1:dim(p)[2]){
  lines(y=p[,i], x=seq(-ra,ra,length.out=(ra*2+1)*10), ...)
}
#pos <- (1:((ra*2+1)*10))[((1:length(-ra:ra))*10)-ra]

ticks <- unique(round(seq(-ra,ra,length.out=(ra*2+1)*10)))

length(ticks)
length(-ra:ra)
axis(1,ticks,labels=as.character(-ra:ra))

abline(h = 0, v = na.omit(pair_obj$threshold[itemnumber,]), col = "gray60")

hohe <- apply(p,2,max)
lage <- apply(p,2,function(x){as.numeric(rownames(p)[max(x)==x])  }    )

text(lage,hohe,labels=names(lage),pos=1,cex=.7)
}
##### silent return -------------------------

if (plot==FALSE){(return(p))}
}


