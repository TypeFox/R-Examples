#' Expected Maximum Difference
#' @param x a column vector of scores on which the rsd is conditioned
#' @param o a column vector of equated scores based on the overall population (aligned with elements in x)
#' @param g column vectors of equated scores based on various subpopulations (aligned with elements in x)
#' @param n a scalar indicating the number of groups
#' @param f a column vector of relative frequency associated with each raw score (can be based on either overall population or a subpopulation) (aligned with elements in x)
#' @param s a scalar representing the standard deviation of x for any (sub)population of interest (e.g., synthetic population) (default is 1, which leads to calculation of the unstandardized emaxd)
#' @author Anne Corinne Huggins-Manley 
#' @references 
#' \itemize{
#' \item{Dorans, N.J., & Holland, P.W. (2000). Population invariance and the equitability of tests: Theory and the linear case. Journal of Educational Measurement, 37, 281-306.}
#' }
#' @seealso \code{\link{maxd}}
#' @examples
#' #Unstandardized EMAXD across subpopulation 1 and subpopulation 2 in the example data set, ex.data
#' emaxd(x=ex.data[,1],o=ex.data[,2],g=c(ex.data[,3],ex.data[,4]),n=2,f=ex.data[,8])
#' 
#' #Unstandardized EMAXD across subpopulations 1 thru 5 in the example data set, ex.data
#' emaxd(x=ex.data[,1],o=ex.data[,2],g=c(ex.data[,3],ex.data[,4],ex.data[,5],ex.data[,6],ex.data[,7]),n=5,f=ex.data[,8])
#' 
#' #Standardized EMAXD across subpopulations 1 thru 5 in the example data set, ex.data
#' emaxd(x=ex.data[,1],o=ex.data[,2],g=c(ex.data[,3],ex.data[,4],ex.data[,5],ex.data[,6],ex.data[,7]),n=5,f=ex.data[,8],s=4.2)
#' @return expected maximum difference
#' @export

emaxd<- function(x,o,g,n,f,s){
  
  if(missing(s))
    s <- 1
  
  dm <- c() 
  t<-as.matrix(g)
  dim(t)<-c(length(x),n)  
  
  for (i in 1:n){
    d <- assign(x=paste("rsd",i,sep=""),value=sqrt((o-t[,i])^2))
    
    dm <- cbind(dm,d)  
  }
  
  mx<-c()
  
  for (i in 1:length(x)){
    m<- assign(x=paste("max",i,sep=""),value=max(dm[i,]))
    mx <- rbind(mx,m,row.names=NULL)
  }
  
  mx2 <- mx*f
  mx3 <- (sum(mx2))/s
  
  return(mx3)
}