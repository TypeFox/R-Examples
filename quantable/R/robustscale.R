#' robust scaling
#' uses median an mad instead of mean and row
#' applies the scaling to the columns (samples) by default
#' @export
#' @param data matrix or data.frame
#' @param dim - should rows (1) or columns (2:default) be scaled
#' @param center - subract median (default:TRUE)
#' @param scale - scale by mad  (default:FALSE)
#' @examples
#' tmp = matrix(rep((1:100),times = 4) + rnorm(100*4,0,3),ncol=4)
#' mean = c(20,30,10,40)
#' sd = c(4,3,4,5)
#' tmp = sweep(tmp,2,sd,"*")
#' tmp = sweep(tmp,2,mean,"+")
#' boxplot(tmp)
#' tmp = robustscale(tmp)
#' boxplot(tmp$data)
robustscale <- function(data, dim=2, center=TRUE, scale=TRUE){
  medians = NULL
  if(center){
    medians <- apply(data,dim,median,na.rm=TRUE)
    data = sweep(data,dim,medians,"-")
  }
  mads=NULL
  if(scale){
    mads <- apply(data,dim, mad,na.rm =TRUE)
    data = (sweep(data,dim,mads,"/"))
  }
  return(list(data=data,medians=medians,mads=mads))
}
