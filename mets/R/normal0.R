
##' @export
loglikMVN <- function(yl,yu,status,mu,S,thres) {
    .Call("loglikMVN",
          yl=as.matrix(yl),
          yu=as.matrix(yu),
          status=as.integer(status),
          mu=as.matrix(mu),dmu=NULL,s=as.matrix(S),ds=NULL,
          z=NULL,su=NULL,dsu=NULL,
          threshold=as.matrix(thres),dthreshold=NULL,
          score=FALSE,
          package="mets")        
}


##' @export
scoreMVN <- function(y,mu,S,dmu,dS) {
    .Call("loglikMVN",
          yl=as.matrix(y),
          yu=NULL,
          status=NULL,
          mu=as.matrix(mu),dmu=as.matrix(dmu),
          s=as.matrix(S),ds=as.matrix(dS),
          z=NULL,su=NULL,dsu=NULL,
          threshold=NULL,dthreshold=NULL,
          score=TRUE,
          package="mets")        
}
