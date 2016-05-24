#' Mean rotation
#'
#' Compute the sample geometric or projected mean.
#'
#' This function takes a sample of 3D rotations (in matrix or quaternion form) and returns the projected arithmetic mean denoted \eqn{\widehat{\bm S}_P}{S_P} or
#' geometric mean \eqn{\widehat{\bm S}_G}{S_G} according to the \code{type} option.
#' For a sample of \eqn{n} rotations in matrix form \eqn{\bm{R}_i\in SO(3), i=1,2,\dots,n}{Ri in SO(3), i=1,2,\dots,n}, the mean-type estimator is defined as \deqn{\widehat{\bm{S}}=argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd^2(\bm{R}_i,\bm{S})}{argmin\sum d^2(Ri,S)} 
#' where \eqn{d} is the Riemannian or Euclidean distance.
#' For more on the projected mean see \cite{moakher02} and for the geometric mean see \cite{manton04}.
#' For the projected mean from a quaternion point of view see \cite{tyler1981}.
#' 
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix form (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param type string indicating "projected" or "geometric" type mean estimator.
#' @param epsilon stopping rule for the geometric-mean.
#' @param maxIter maximum number of iterations allowed for geometric-mean.
#' @param ... additional arguments.
#' @return Estimate of the projected or geometric mean of the sample in the same parameterization.
#' @aliases mean.Q4
#' @seealso \code{\link{median.SO3}}, \code{\link{bayes.mean}}, \code{\link{weighted.mean.SO3}}
#' @@cite tyler1981, moakher02, manton04
#' @export 
#' @method mean SO3
#' @examples
#' Rs <- ruars(20, rvmises, kappa = 0.01)
#' mean(Rs)                               #Projected mean
#' project.SO3(colMeans(Rs))              #Same as mean(Rs)
#' 
#' mean(Rs, type = "geometric")           #Geometric mean
#' rot.dist(mean(Rs))                     #Bias of the projected mean
#' rot.dist(mean(Rs, type = "geometric")) #Bias of the geometric mean
#' 
#' Qs <- as.Q4(Rs)
#' mean(Qs)                               #Projected mean
#' mean(Qs, type = "geometric")           #Geometric mean
#' rot.dist(mean(Qs))                     #Bias of the projected mean
#' rot.dist(mean(Qs, type = "geometric")) #Bias of the geometric mean

mean.SO3 <- function(x, type = "projected", epsilon = 1e-05, maxIter = 2000, ...) {
	
	Rs<-formatSO3(x)	
	
	if(nrow(Rs)==1)
		return(Rs)
	
  type <- try(match.arg(type,c("projected", "geometric")),silent=T)
  
  if (class(type)=="try-error")
    stop("type needs to be one of 'projected' or 'geometric'.")

	
	if(type=='projected'){
		R<-meanSO3C(Rs)
	}else{
		R<-gmeanSO3C(Rs,maxIter,epsilon)
	}
	
  class(R)<-"SO3"
  return(R)
}

#' @rdname mean.SO3
#' @aliases mean.SO3
#' @export 
#' @method mean Q4

mean.Q4 <- function(x, type = "projected", epsilon = 1e-05, maxIter = 2000,...) {
	
	Qs<-formatQ4(x)
	
	if(nrow(Qs)==1)
		return(Qs)
	
	type <- try(match.arg(type,c("projected", "geometric")),silent=T)
	
	if (class(type)=="try-error")
	  stop("type needs to be one of 'projected' or 'geometric'.")
  
	if(type=='projected'){
		
		R<-meanQ4C(Qs)
		
	}else{
		
		Rs<-as.SO3.Q4(Qs)
  	R<-gmeanSO3C(Rs,maxIter,epsilon)
		R<-as.Q4.SO3(R)
	}
	
	class(R)<-'Q4'
  return(R)
  
}


#' Median rotation
#' 
#' Compute the sample projected or geometric median.
#'
#' The median-type estimators are defined as \deqn{\widetilde{\bm{S}}=argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd(\bm{R}_i,\bm{S}).}{argmin\sum d(Ri,S).}  
#' If the choice of distance metric \eqn{d} is Riemannian then the estimator is called the geometric median, 
#' and if the distance metric in Euclidean then it is called the projected median.
#' The algorithm used in the geometric case is discussed in \cite{hartley11} 
#' and the projected case is in \cite{stanfill2013}.
#'
#' @name median.SO3
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix form (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param type string indicating "projected" or "geometric" type mean estimator.
#' @param epsilon stopping rule.
#' @param maxIter maximum number of iterations allowed before returning most recent estimate.
#' @param ... additional arguments.
#' @return Estimate of the projected or geometric median in the same parameterization.
#' @aliases median.Q4 median.SO3
#' @seealso \code{\link{mean.SO3}}, \code{\link{bayes.mean}}, \code{\link{weighted.mean.SO3}}
#' @@cite hartley11 stanfill2013
#' @export
#' @examples
#' Rs <- ruars(20, rvmises, kappa = 0.01)
#' median(Rs)                               #Projected median
#' median(Rs, type = "geometric")           #Geometric median
#' rot.dist(median(Rs))                     #Bias of the projected median
#' rot.dist(median(Rs, type = "geometric")) #Bias of the geometric median
#' 
#' Qs <- as.Q4(Rs)
#' median(Qs)                               #Projected median
#' median(Qs, type = "geometric")           #Geometric median
#' rot.dist(median(Qs))                     #Bias of the projected median
#' rot.dist(median(Qs, type = "geometric")) #Bias of the geometric median

median<-function(x,...){
  UseMethod("median")
}

#' @rdname median.SO3
#' @aliases median.Q4 median
#' @method median SO3
#' @export 

median.SO3 <- function(x, type = "projected", epsilon = 1e-05, maxIter = 2000,...) {
  
	Rs<-formatSO3(x)
	n<-nrow(Rs)
	
	if(nrow(Rs)==1)
		return(Rs)
	
	type <- try(match.arg(type,c("projected", "geometric")),silent=T)
	
	if (class(type)=="try-error")
	  stop("type needs to be one of 'projected' or 'geometric'.")
  
  if (type == "projected") {
		
  	S<-medianSO3C(Rs,maxIter,epsilon)
      
  } else {
      
		S<-HartmedianSO3C(Rs,maxIter,epsilon)
      
  }
    
	class(S)<-"SO3"
  return(S)
}


#' @rdname median.SO3
#' @aliases median.SO3 median
#' @method median Q4
#' @export 

median.Q4 <- function(x, type = "projected", epsilon = 1e-05, maxIter = 2000,...) {
	
	Qs<-formatQ4(x)
	
	if(length(Qs)==4)
		return(Qs)

  Rs<-as.SO3.Q4(Qs)
  
  R<-median.SO3(Rs,type,epsilon,maxIter,...)

  return(as.Q4.SO3(R))
}


#' Weighted mean rotation
#'
#' Compute the weighted geometric or projected mean of a sample of rotations.
#'
#' This function takes a sample of 3D rotations (in matrix or quaternion form) and returns the weighted projected arithmetic mean \eqn{\widehat{\bm S}_P}{S_P} or
#' geometric mean \eqn{\widehat{\bm S}_G}{S_G} according to the \code{type} option.
#' For a sample of \eqn{n} rotations in matrix form \eqn{\bm{R}_i\in SO(3), i=1,2,\dots,n}{Ri in SO(3), i=1,2,\dots,n}, the weighted mean is defined as 
#' \deqn{\widehat{\bm{S}}=argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nw_id^2(\bm{R}_i,\bm{S})}{argmin\sum wi d^2(Ri,S)} where \eqn{d}
#' is the Riemannian or Euclidean distance.  For more on the projected mean see \cite{moakher02} and for the geometric mean see \cite{manton04}.
#'
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix form (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param w vector of weights the same length as the number of rows in x giving the weights to use for elements of x.
#' @param type string indicating "projected" or "geometric" type mean estimator.
#' @param epsilon stopping rule for the geometric method.
#' @param maxIter maximum number of iterations allowed before returning most recent estimate.
#' @param ... only used for consistency with mean.default.
#' @return Weighted mean of the sample in the same parameterization.
#' @seealso \code{\link{median.SO3}}, \code{\link{mean.SO3}}, \code{\link{bayes.mean}}
#' @aliases weighted.mean.Q4
#' @@cite moakher02
#' @export 
#' @method weighted.mean SO3
#' @examples
#' Rs <- ruars(20, rvmises, kappa = 0.01)
#' mean(Rs)                   #Find the equal-weight projected mean
#' wt <- abs(1/mis.angle(Rs)) #Use the rotation misorientation angle as weight
#' weighted.mean(Rs, wt)      #as weight
#' rot.dist(mean(Rs))
#' rot.dist(weighted.mean(Rs, wt)) #usually much smaller than unweighted mean
#' 
#' #Can do the same thing with quaternions
#' Qs <- as.Q4(Rs)
#' mean(Qs) 
#' wt <- abs(1/mis.angle(Qs)) 
#' weighted.mean(Qs, wt)      
#' rot.dist(mean(Qs))
#' rot.dist(weighted.mean(Qs, wt)) 

weighted.mean.SO3 <- function(x, w, type = "projected", epsilon = 1e-05, maxIter = 2000, ...) {
	
	Rs<-formatSO3(x)
	
	if(nrow(Rs)==1)
		return(Rs)
	
	if(length(w)!=nrow(Rs))
		stop("'Rs' and 'w' must have same length")
	
	type <- try(match.arg(type,c("projected", "geometric")),silent=T)
	
	if (class(type)=="try-error")
	  stop("type needs to be one of 'projected' or 'geometric'.")
	
	if(any(w<0))
		warning("Negative weights were given.  Their absolute value is used.")
	
	w<-abs(w/sum(w))
	
	wRs<-w*Rs
	
	R <- as.SO3(project.SO3(matrix(colSums(wRs), 3, 3)))
	
	if (type == "geometric") {
		n <- nrow(Rs)
		d <- 1
		iter <- 0
		s <- matrix(0, 3, 3)
		
		while (d >= epsilon) {
			
			R <- R %*% skew.exp(s)
			
			s <- matrix(colSums(w*t(apply(Rs, 1, tLogMat, S = R))), 3, 3)
			
			d <- norm(s, type = "F")
			
			iter <- iter + 1
			
			if (iter >= maxIter) {
				warning(paste("No convergence in ", iter, " iterations."))
        class(R)<-"SO3"
				return(R)
			}
		}
    class(R)<-"SO3"
		#R<-as.SO3(R)	
	}
	
	return(R)
}

#' @rdname weighted.mean.SO3
#' @aliases weighted.mean.SO3
#' @export 
#' @method weighted.mean Q4


weighted.mean.Q4 <- function(x, w, type = "projected", epsilon = 1e-05, maxIter = 2000,...) {
	
	Qs<-formatQ4(x)
	
	if(nrow(Qs)==1)
		return(Qs)
	
	Rs<-as.SO3(Qs)
	
	R<-weighted.mean.SO3(Rs,w,type,epsilon,maxIter)
	
	return(as.Q4.SO3(R))
	
}
