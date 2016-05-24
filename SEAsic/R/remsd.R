#' Root Expected Mean Square Difference
#' @param x a column vector of scores on which the rmsd is conditioned
#' @param o a column vector of equated scores based on the overall population (aligned with elements in x)
#' @param g column vectors of equated scores based on various subpopulations (aligned with elements in x)
#' @param f a column vector of relative frequency associated with each raw score (can be based on either overall population or a subpopulation) (aligned with elements in x)
#' @param s a scalar representing the standard deviation of x for any (sub)population of interest (e.g., synthetic population) (default is 1, which leads to calculation of the unstandardized remsd)
#' @param w A row vector of weights for subpopulations 1 thru n (length = number of groups)
#' @author Anne Corinne Huggins-Manley 
#' @references 
#' \itemize{
#' \item{Dorans, N.J., & Holland, P.W. (2000). Population invariance and the equitability of tests: Theory and the linear case. Journal of Educational Measurement, 37, 281-306.}
#' }
#' @seealso \code{\link{rmsd}}
#' @examples
#' #Unstandardized REMSD for subpopulations 1 and 2 in the example data set, ex.data, assuming equal weights for the subpopulations
#' remsd(x=ex.data[,1],o=ex.data[,2],g=c(ex.data[,3],ex.data[,4]),f=ex.data[,8],w=c(.5,.5))
#' 
#' #Unstandardized REMSD for all five subpopulations in the example data set, ex.data
#' remsd(x=ex.data[,1],o=ex.data[,2],g=c(ex.data[,3],ex.data[,4],ex.data[,5],ex.data[,6],ex.data[,7]),f=ex.data[,8],w=c(.1,.2,.4,.2,.1))
#' 
#' #Standardized REMSD for all five subpopulations in the example data set, ex.data
#' remsd(x=ex.data[,1],o=ex.data[,2],g=c(ex.data[,3],ex.data[,4],ex.data[,5],ex.data[,6],ex.data[,7]),f=ex.data[,8],w=c(.1,.2,.4,.2,.1),s=4.2)
#' @note The equally weighted version of this index (Kolen & Brennan, 2004) can be obtained by inputting a w vector consisting of identical elements that sum to 1. See example 1 above.
#' @return root expected mean square difference
#' @export

remsd<- function(x,o,g,f,s,w){
  tm <- c() 
  
  for (i in 1:length(w)){
    t<- assign(x=paste("t",i,sep=""),value=c(rep(w[i],length(x))))
    
    tm <- c(tm,t,row.names=NULL)
  }
  
  if(missing(s))
    s <- 1
  
  rmsdX<-(tm*(o-g)^2)
  rmsdY<-rowSums(matrix(rmsdX,nrow=length(x))) 
  remsdX<-rmsdY*f
  remsdY<-sum(remsdX)
  remsdZ<-(sqrt(remsdY))/s  
  return(remsdZ)
}