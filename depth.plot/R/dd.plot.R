#' @title Depth-Depth Plots
#'
#' @description \code{dd.plot} is a multivariate genralization of a normal \code{QQ-plot}.
#'    It produces a DD-plot of two datasets.
#' @param data1 A matrix or a data.frame with each row as a p-variate observation.
#' @param data2 A matrix or a data.frame (defaults to a standard independent p-variate normal).
#' @param main Plot labels. The title of the plot.
#' @param xlab Plot labels. The \code{x-axis} label of the plot.
#' @param ylab Plot labels. The \code{y-axis} label of the plot.
#' @param col The color of the points
#' @param pch character string or vector of 1-characters or integers for plotting characters.
#' @import mvtnorm
#' @import graphics
#' @author Somedip Karmakar <somedip@yahoo.co.in>
#' @author Omker Mahalanobish <omker.scorpio@gmail.com>
#' @export
#' @seealso \code{\link{spatial.depth}}
#' @return A \code{DD-plot} of the input data
#' @examples
#' u<-matrix(rnorm(300,1,4),ncol=3)
#' dd.plot(u)
#'
dd.plot<-function(data1,data2=rmvnorm(nrow(data1),array(0,ncol(data1)),diag(1,ncol(data1),ncol(data1))),
                  main="Normal DD-plot",xlab="Sample Depths",ylab="Normal Depths",col="black",pch=20){
  spatial.depth=function(data,x){
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
    sd
  }
  something1<-c(spatial.depth(data1,data1),spatial.depth(data1,data2))
  something2<-c(spatial.depth(data2,data1),spatial.depth(data2,data2))
  x<-cbind(something1,something2)
  matplot(x[,1],x[,2],main=main,xlab=xlab,ylab=ylab,col=col,pch=pch)
}
