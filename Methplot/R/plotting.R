#' This function visualizes the dataframe yielded by function "getdata".
#'
#' Prerequisites: You need to install packages: "ggplot2", "grid", and "reshape"
#'
#' @param x the dataframe name that is in the format of "getdata" output
#' @param x.title the title of the output plot
#' @param condition a character vector that gives the sample list that you want to plot. These samples will be merged together by adding up the read numbers to yield one plot.
#' @param n the number of CpG sites in this region.
#' @param legendpos the position of the legend. Default: "null" ("left", "right", "bottom", "top", or two-element numeric vector) 
#' @return Output plot: Rows indicate reads and are presented in percentage (y axis), showing whether each of CpG sites (x axis) in the target region is methylated (light green) or demethylated (dark green). Reads are sorted so that those with no demethylated positions are at the bottom, and those with most demethylated positions are at the top.
#' @export
#' @author Xin Yang \email{xin.yang@@cimr.cam.ac.uk}
#' @examples
#' data(mydata)
#' plotdata(mydata, x.title="Methylation Plot", condition="P1_A1", n=10, legendpos="right")
#' @import ggplot2 grid reshape
if(getRversion() >= "2.15.1")  utils::globalVariables(c("x1", "x2","y1","y2","base"))

plotdata<- function(x,x.title, condition,n, legendpos="null") {
  	
  if(length(condition)>1) {
    x<-apply(x[,condition], 1, sum)
  }else{
    x<-x[,condition]}
  x<-as.data.frame(x)
  x$nT <- nchar(gsub("C","",rownames(x)))

  x$prev <- strReverse(rownames(x))
  x <- x[order(x$nT,x$prev),]
  ss<-strsplit(as.character(rownames(x)),"")
  ss <- t(as.data.frame(ss))
  colnames(ss) <- sprintf("s%s",1:n)
  rownames(ss) <- NULL
  x <- cbind(x,ss)

  x$y2 <- cumsum(x[,1])
  x$y2 <- 100*x$y2/max(x$y2)
  x$y1 <- c(0,x$y2[-nrow(x)])
  df <- lapply(1:n, function(i) {
    data.frame(x1=i-1, x2=i, y1=x$y1, y2=x$y2, base=x[,sprintf("s%s",i)])
  })
  df <- do.call("rbind",df)
  
  head(x2 <- melt(x[,c("nT","y1","y2")], id="nT"))

  ggplot(data=df, aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=base)) +
    geom_rect() +
    scale_fill_brewer(palette="Greens") + xlab("position") +
    ylab(paste("proportion (total reads = ", sum(x[,1]),")", sep="")) +
    labs(title=x.title)+
    theme(legend.position=legendpos)+
    scale_x_discrete(breaks=1:n, labels=sprintf("s%s", 1:n), expand=c(0.05, -0.4))+
    theme_set(theme_bw())

}

strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste,
                                 collapse="")
