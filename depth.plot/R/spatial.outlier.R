#' @title Multivariate Spatial Outlier
#'
#' @description \code{spatial.outlier} is used to find the multivariate spatial outlier within a p-variate data cloud or to identify if any p-variate observation is an outlier with respect to a p-variate data cloud.
#' @param data A matrix or a data.frame of p-variate observations which works as the data cloud.
#' @param x A matrix or a data.framep-variate to test whether is an outlier with respect to the \code{data}. Defaults to \code{data}, to find outliers (if exists) within the data.
#' @param threshold A decimal threshold between \code{0} and \code{1} on the \code{\link{spatial.depth}}. \code{Spatial depth} values less than which will be considered as \code{outlier}. Defaults to \code{0.05}. Usually taken as \code{0.1} or \code{0.05} or \code{0.01}.
#' @author Somedip Karmakar <somedip@yahoo.co.in>
#' @author Omker Mahalanobish <omker.scorpio@gmail.com>
#' @export
#' @return \code{FALSE}  :: If there doesnot exist any outlier
#' @return A list with objects (If outliers exist)
#' @return \code{index}  :: Returns the indices of the outliers
#' @return \code{observation}  :: Returns the p-variate outliers
#' @examples
#' u<-matrix(rnorm(60,0,1),ncol=3)
#' u0<-matrix(runif(9,3,4),ncol=3)
#' spatial.outlier(u,rbind(u,u0))

spatial.outlier<-function(data,x=data,threshold=0.05){
  spatial.depth=function(x,data){
    if(ncol(as.matrix(x))==1)
      x<-t(x)
    if(ncol(data)!=ncol(as.matrix(x))){
      sd<-"Dimensions do not match"
    }else{
      spd<-function(data,x)
      {
        v=array(0,ncol(data))
        for(j in 1:ncol(data))
        {
          for(i in 1:nrow(data))
          {
            if(sqrt(sum((x-data[i,])^2))!=0)
              v[j]=v[j]+((x[j]-data[i,j])/sqrt(sum((x-data[i,])^2)))
          }
          v[j]=v[j]/nrow(data)
        }
        sd=1-sqrt(sum(v^2))
      }
      sd<-apply(x,1,function(y){spd(data,y)})
    }
    sd
  }
  sd<-spatial.depth(x,data)
  if(length(dim(x))!=0){
    m<-list()
    m[[1]]<-which(sd<=threshold)
    m[[2]]<-x[m[[1]],]
    names(m)<-c("index","observation")
    if(length(m[[1]])==0)
      m<-"FALSE"
  }else{
    m<-ifelse(sd<=threshold,"TRUE","FALSE")
  }
  return(m)
}
