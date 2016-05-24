#' @title Spatial Depth
#'
#' @description \code{spatial.depth} is used to find the spatial depth of one or more p-variate observation(s) in a data cloud of numerous p-variate observations.
#' @param x A matrix or a data.frame of objects (numerical vector as one object) whose depth is to be found; each row consists a p-variate observation.
#' @param data A matrix or a data.frame of objects which acts as the data cloud. Each row consists of a p-variate observation.
#' @author Somedip Karmakar <somedip@yahoo.co.in>
#' @author Omker Mahalanobish <omker.scorpio@gmail.com>
#' @export
#' @return Numerical vector of depths, one for each row in \code{x}; or one depth value if \code{x} is numerical.
#' @examples
#' u<-matrix(rnorm(90,0,1),ncol=3)
#' u0<-matrix(runif(9,0,1),ncol=3)
#' spatial.depth(u0,u)

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
  print(sd)
}
