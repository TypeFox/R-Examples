#' @title Spatial Rank
#'
#' @description Used to compute the \code{Spatial Rank} of a p-variate observation with respect to a p-variate data cloud.
#' @param x A numeric p-variate vector whose spatial rank is to be calculated.
#' @param data A matrix or a data.frame with each row as a p-variate observation.
#' @author Somedip Karmakar <somedip@yahoo.co.in>
#' @author Omker Mahalanobish <omker.scorpio@gmail.com>
#' @export
#' @return The spatial rank of \code{x} with respect to \code{data}.
#' @examples
#' u<-matrix(rnorm(90,0,1),ncol=3)
#' u0<-runif(3,0,1)
#' spatial.rank(u0,u)

spatial.rank=function(x,data)
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
  v
}

#' @title Multivariate Quantile
#'
#' @description Used to compute the \code{p-variate quantile} of a p-variate observation with respect to a p-variate data cloud.
#' @param x A numeric p-variate \code{spatial rank}. Elements must lie within -1 and +1, with a 0-vector denoting the median.
#' @param data A matrix or a data.frame with each row as a p-variate observation.
#' @author Somedip Karmakar <somedip@yahoo.co.in>
#' @author Omker Mahalanobish <omker.scorpio@gmail.com>
#' @import stats
#' @export
#' @return The \code{x}th mutivariate quantile with respect to \code{data}.
#' @seealso \code{\link{spatial.rank}}
#' @examples
#' u<-matrix(rnorm(90,0,1),ncol=3)
#' u0<-runif(3,0,1)
#' multi.quant(spatial.rank(u0,u),u)

multi.quant=function(x,data)
{
  fr <- function(y) {
    s=0
    for(i in 1:nrow(data))
    {
      s=s+sqrt(sum((data[i,]-y)^2))+t(x)%*%(data[i,]-y)
    }
    s
  }
  zq=optim(rep(0,ncol(data)), fr)$par
  zq
}


