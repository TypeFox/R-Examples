#' Median polish multidimensional.
#'
#' Fits an additive model for multidimensional array, using Tukey's median polish procedure.
#' 

#' @param data object of class array, table, or matrix, see details.
#' @param eps real number greater than \code{0}, default 0.01. A tolerance for convergence: see Details
#' @param maxiter the maximum number of iterations. Default 10.
#' @param na.rm logical. If the data contains NA's. Default TRUE.
#' @param \dots ignored.
#' @details the model fitted is additive \eqn{constant + dim_{1} + dim_{2} + \cdots + dim_{n}}. The algorithm works by alternately removing medians of \eqn{dim_{1}, \cdots, dim_{n}}, and continues until the proportional reduction in the sum of absolute residuals is less than eps or until there have been maxiter iterations. If na.rm is FALSE the presence of any NA value in x will cause an error, otherwise NA values are ignored. MedianPolishM returns an object of class MedianPolishM (see below). There are a plotting method for this class, \code{\link{plot.MedianPolishM}}.
#' @return An object of class medpolish with the following named components in a list:
#' @return \item{residuals}{the residuals.}
#' @return \item{overall}{the fitted constant term.}
#' @return \item{effects}{the fitted every dimensions effects of array multidimensional.}
#' @return \item{iter}{number of iterations used in the range maxiter.}
#' @references Hoaglin, D. C., Mosteller, F., & Tukey, J. W. (Eds.). (2011). \emph{Exploring data tables, trends, and shapes} (Vol. 101). John Wiley & Sons.\href{http://www.wiley.com/WileyCDA/WileyTitle/productCd-047004005X.html}{[link]}
#' @examples A<-MedianPolishM(UCBAdmissions, eps=0.1, maxiter=2, na.rm=TRUE)
#' plot(A)
#' @export 
#' 
MedianPolishM.default <-
  function(data, eps=0.01, maxiter=10L, na.rm=TRUE,...){   
    
    stopifnot(inherits(data, c("array","table","matrix")))
    data0<-data
    Res<-function(x,w){ aperm( aperm(x,c(w,c(1:length(dim(x)))[-w]))-apply(x,w,median,na.rm=na.rm), order(c(w,c(1:length(dim(x)))[-w])) ) }
    
    Effect<-list()
    Effect01<-list()
    Effect1<-list()
    Effect2<-list()
    Overall=0
    
    for(qu in 1L:maxiter){
      
      for(i in 1:length(dim(data))){
        Effect[[i]]<-apply(data,i,median,na.rm=na.rm)
        Effect01[[i]]<-rep(0,dim(data)[i])
        data<-Res(data,i)
      }
      if(qu==1)
      {
        Effect1<-Effect01
        for(i in 1:length(dim(data))){
          Effect1[[i]]<-Effect1[][[i]]+Effect[][[i]]}
      }else
      {
        for(i in 1:length(dim(data))){
          Effect1[[i]]<-Effect1[][[i]]+Effect[][[i]]}
      }
      
      converged<-matrix(0,1,(length(dim(data))))
      for(i in 1:length(dim(data))){
        converged[i]<-abs(median(apply(data, i,median,na.rm=na.rm),na.rm=na.rm))<= eps
      }
      
      if (prod(converged)==1) 
        break
    } 
    
    for(j in 1:length(dim(data))){
      Effect2[[j]]<-Effect1[[j]]-median(Effect1[[j]],na.rm=na.rm)
      Overall<-median(Effect1[[j]],na.rm=na.rm)+Overall}
    names(Effect2)<-names(dimnames(data))
    
    MedPolish <- list(residuals=data,overall=Overall,effects=Effect2,iter=qu,data0=data0) 
    MedPolish$data<-data0
    MedPolish$eps<-eps
    MedPolish$Gr<-1
    
    class(MedPolish) <- "MedianPolishM"
    return(MedPolish)
    
  }

