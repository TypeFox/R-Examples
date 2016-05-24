#' Maximum Difference
#' @param x a column vector of scores on which the maximum difference is conditioned
#' @param o a column vector of equated scores based on the overall population (aligned with elements in x)
#' @param g column vectors of equated scores based on various subpopulations (aligned with elements in x)
#' @param n a scalar indicating the number of groups
#' @param d a scalar of the difference that matters
#' @param s a scalar representing the standard deviation of x for any (sub)population of interest (e.g., synthetic population) (default is 1, which leads to calculation of the unstandardized maxd)
#' @param ymax A maximum value for the y axis of the plot (default = 3 + the maximum maxD value)
#' @param xlab A label for the x axis of the plot (default = Score Scale)
#' @param color color of maxD line and points on plot (default = red)
#' @author Anne Corinne Huggins-Manley 
#' @references 
#' \itemize{
#' \item{Dorans, N.J., & Holland, P.W. (2000). Population invariance and the equitability of tests: Theory and the linear case. Journal of Educational Measurement, 37, 281-306.}
#' }
#' @seealso \code{\link{emaxd}}
#' @examples
#' #Unstandardized MAXD across subpopulation 1 and subpopulation 2 in the example data set, ex.data
#' maxd(x=ex.data[,1],o=ex.data[,2],g=c(ex.data[,3],ex.data[,4]),n=2,d=.5)
#' 
#' #Unstandardized MAXD across subpopulations 1 thru 5 in the example data set, ex.data
#' maxd(x=ex.data[,1],o=ex.data[,2],g=c(ex.data[,3],ex.data[,4],ex.data[,5],ex.data[,6],ex.data[,7]),n=5,d=.5)
#' 
#' #Standardized MAXD across subpopulations 1 thru 5 in the example data set, ex.data
#' maxd(x=ex.data[,1],o=ex.data[,2],g=c(ex.data[,3],ex.data[,4],ex.data[,5],ex.data[,6],ex.data[,7]),n=5,d=.5,s=4.2)
#' @return A data frame of maximum difference indices, conditioned on the score scale
#' @return A plot of the maximum difference indices in reference to the difference that matters
#' @export

maxd<- function(x,o,g,n,d,s,ymax,xlab,color){
  
  if(missing(s))
    s <- 1
  
  dm <- c() 
  t<-as.matrix(g)
  dim(t)<-c(length(x),n)  
  
  for (i in 1:n){
    dn <- assign(x=paste("rsd",i,sep=""),value=sqrt((o-t[,i])^2))
    
    dm <- cbind(dm,dn)  
  }
  
  mx<-c()
  
  for (i in 1:length(x)){
    m<- assign(x=paste("max",i,sep=""),value=max(dm[i,]))
    mx <- rbind(mx,m,row.names=NULL)
  }
  mx2 <- cbind(x,mx/s)
  
  mx3<- data.frame(mx2,row.names=NULL)
  colnames(mx3)<- c("x","maxd")
  
  if(missing(ymax))
    ymax <- max((mx3[,2])+3)
  if(missing(xlab))
    xlab <- "Score Scale"
  if(missing(color))
    color <- "red"
  
  
  plot<- plot(type="b",x=mx3[,1],y=mx3[,2],ylim=c(0,ymax),xlab=xlab,ylab="Maximum Difference",col=color,pch=19)
  abline(d/s,0,col="black",lty=2)
  legend(x=min(mx3[,1]),y=ymax,c("MaxD","DTM"),col=c(color,"black"),pch=c(19,-1),lty=c(0,5))
  
  return(mx3)
  return(plot)
}