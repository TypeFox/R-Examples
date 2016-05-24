#' K-means algorithm for the clustering of variables
#'
#' K-means algorithm for the clustering of variables. Directional or local groups may be defined. 
#' Each group of variables is associated with a latent component. 
#' Moreover external information collected on the observations or on the variables may be introduced.
#' 
#' The initalization can be made at random, repetitively, or can be defined by the user.
#' 
#' The parameter "strategy" makes it possible to choose a strategy for setting aside variables
#' that do not fit into the pattern of any cluster.   
#'
#' @param X The matrix of the variables to be clustered
#' @param Xu The external variables associated with the columns of X
#' @param Xr The external variables associated with the rows of X
#' @param method The criterion to use in the cluster analysis.\cr 
#'        1 or "directional" : the squared covariance is used as a measure of proximity (directional groups). \cr    
#'        2 or "local"       : the covariance is used as a measure of proximity (local groups)
#' @param sX TRUE/FALSE : standardization or not of the columns X (TRUE by default)\cr
#'        (predefined -> cX = TRUE : column-centering of X)
#' @param sXr TRUE/FALSE : standardization or not of the columns Xr (FALSE by default)\cr
#'        (predefined -> cXr    = TRUE : column-centering of Xr)
#' @param sXu TRUE/FALSE : standardization or not of the columns Xu (FALSE by default)\cr
#'        (predefined -> cXu= FALSE : no centering, Xu considered as a weight matrix)
#' @param clust : a number i.e.  the size of the partition, K,
#'        or  a vector of INTEGERS i.e. the group membership of each variable in the initial partition (integer between 1 and K)
#' @param iter.max maximal number of iteration for the consolidation (20 by default)
#' @param nstart nb of random initialisations in the case where init is a number  (100 by default)
#' @param strategy "none" (by default), or "kplusone" (an additional cluster for the noise variables),
#'        or "sparselv" (zero loadings for the noise variables)
#' @param rho a threshold of correlation between 0 and 1 (0.3 by default)
#' 
#' @return \item{tabres}{ 
#'         The value of the clustering criterion at convergence.\cr
#'         The percentage of the explained initial criterion value.\cr
#'         The number of iterations in the partitioning algorithm.}
#'         \item{clusters}{ the group's membership}
#'         \item{comp}{ The latent components of the clusters}
#'         \item{loading}{ if there are external variables Xr or Xu :  The loadings of the external variables}
#' @seealso CLV, LCLV
#' @examples data(apples_sh)
#' #local groups with external variables Xr 
#' resclvkmYX <- CLV_kmeans(X = apples_sh$pref, Xr = apples_sh$senso,method = "local",
#'           sX = FALSE, sXr = TRUE, clust = 2, nstart = 20)
#' @export                        
#'                          
CLV_kmeans <- function(X,Xu=NULL,Xr=NULL,method,sX=TRUE,sXr=FALSE,sXu=FALSE,
                       clust, iter.max=20, nstart=100,strategy="none",rho=0.3)
{
  
  if (method=="directional") method=1
  if (method=="local") method=2
  if(method!=1 & method!=2) stop("method should be 1/directional or 2/local")

  if (missing(clust))
    stop("'clust' must be an integer (the number of clusters) or a vector of vector a vector of intergers (the initial partition)")
  cX=TRUE
  cXr=TRUE
  cXu=FALSE
  
  # verification if some variables have constant values (standard deviation=0)
  who<-which(apply(X,2,sd)==0)
  if ((length(who)>0)&(sX==TRUE)) {
    listwho<-c(": ")
    for (r in 1:length(who)) {listwho=paste(listwho,colnames(X)[who[r]],",")}
    stop("The variables",listwho," have constant values (standard deviation=0). Please remove these variables from the X matrix.")
  }
  if (length(who)>0) {
    listwho<-c(": ")
    for (r in 1:length(who)) {listwho=paste(listwho,colnames(X)[who[r]],",")}
    warning("The variables",listwho," have constant values (standard deviation=0). Please remove these variables from the X matrix.")
  }
  
  X<- scale(X, center=cX, scale=sX)
  p <- ncol(X)
  n <- nrow(X)  

  if (is.null(Xr)) {
    EXTr<-0
  }
  else {
    EXTr<-1                                    
    Xr<- scale(Xr, center=cXr, scale=sXr)
    ntilde <- dim(Xr)[1]
    q<-dim(Xr)[2] 
    if (n != ntilde) stop("X and Xr must have the same number of observations")                              
  }

  if (is.null(Xu)) {
    EXTu<-0
  } 
  else {
    EXTu<-1                   
    Xu<- scale(Xu, center=cXu, scale=sXu) 
    ptilde <- dim(Xu)[1]  
    m<-dim(Xu)[2] 
    if (p != ptilde) {stop("X and Xu must be defined for the same number of
                           variables") }
    if (EXTr==1) {stop("this procedure doesn't allow Xr and Xu to be defined
                       simultenaously. Use l-clv instead")}
  }  
  
  
  crit<-crit_init(method,X,EXTr,Xr,EXTu,Xu)
  sbegin <- sum(crit)  
  
 
 if (length(clust) == 1) {
     K <- clust
     out<-mat_init(X,EXTr,Xr,EXTu,Xu,K)
     comp<-out$comp
     comp <- as.matrix(X[,sort(sample.int(p, K))]) # K columns of X chosen at random
     if (EXTr==1)  {  a<-out$a }
     if (EXTu==1)  {  u<-out$u }
     groupes <- as.factor(consol_affect(method,X,Xr,Xu,EXTr,EXTu,comp,a,u))
  } else {
     nstart = 1
     if (!is.numeric(clust))
         stop("clust must be a vector of integers")
     groupes <- as.factor(clust)
     K <- length(levels(groupes))
     if (p < K)
         stop("more cluster's centers than variables")
     if (length(which(clust > K)) > 0)
         stop("clusters must be numbered from 1 to K")
     if (p != length(groupes))
         stop("the length of clust must be equal to the number of variables")
     out<-mat_init(X,EXTr,Xr,EXTu,Xu,K)
     comp<-out$comp
     if (EXTr==1)  {a<-out$a}
     if (EXTu==1)  {u<-out$u}
 }
           
  if (length((intersect(strategy,c("kplusone","sparselv","none"))))==0) 
      stop("strategy must be either 'kplusone', 'sparselv', 'none'")
 
  if (strategy=="kplusone" & (EXTu!=0 | EXTr!=0) )
       stop(" 'k+1' strategy is not available with external variables, yet")
  if (strategy=="sparselv" & (EXTu!=0 | EXTr!=0) )
       stop(" 'Sparse LV' strategy is not available with external variables, yet")
 
  #####################################################################
  if (strategy=="none" | rho==0 ) 
    listcc = clvk_none(X,n,p,sbegin,EXTr,EXTu,Xu,Xr,method,K,comp,groupes,a,u,iter.max,nstart)
  if (strategy=="sparselv")       
    listcc = clvk_sparse(X,n,p,sbegin,EXTr,EXTu,Xu,Xr,method,K,comp,groupes,a,u,iter.max,nstart,rho)
  if (strategy=="kplusone")       
    listcc = clvk_kp1(X,n,p,sbegin,EXTr,EXTu,Xu,Xr,method,K,comp,groupes,a,u,iter.max,nstart,rho)
  #####################################################################


param<-list(X=X,method=method,n = n, p = p,K = K,nstart = nstart,EXTu=EXTu,EXTr=EXTr,
            sX=sX,sXr=sXr,cXu=cXu,sXu=sXu,strategy=strategy,rho=rho)
listcc= c(listcc, list(param=param)) 


class(listcc) = "clv"
# if (strategy=="sparselv") class(listcc) = "sparselv"
# if (strategy=="kplusone") class(listcc) = "kplusone"

return(listcc)
}
