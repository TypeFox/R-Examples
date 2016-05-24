#' Compute Filter-based Metrics in a Functional Space Between Marked Point Processes
#'
#'
#' The most commonly used and intensively studied metrics for spike trains, which is based on the continuation of event sequence to a real valued continuous function using a smoother function.
#'
#'
#' \code{fmetric} computes filter-based measure between MPP realizations. Discrete event timings are transformed into a continuous function by using a kernel smoother, and usual l2 inner product is adopted for defining the similarity between two point process realizations.
#'
#' 
#' @param S1 marked point process data.
#' @param S2 marked point process data.
#' @param measure \code{"sim"} for similarity and \code{"dist"} for distance. Default \code{"sim"}.
#' @param h filtering function. Default \code{"laplacian"} offers significant computational advantage. A function can be specified here like
#' \code{h=function(x,tau) exp(-x^2/tau)}.
#' The function should be square integrable and non-negative (not checked in the code).
#' @param tau parameter for filtering function.
#' @param M a precision matrix for filter of marks, i.e., exp( - r' M r) is used for filtering marks. It should be symmetric and positive semi-definite.
#' @param abs.tol absolute tolerance for numerical integration.
#' @return Similarity or distance between two inputs (marked) point process S1 and S2.
#' @author Hideitsu Hino \email{hinohide@@cs.tsukuba.ac.jp}, Ken Takano, Yuki Yoshikawa, and Noboru Murata
#' @references M. C. W. van Rossum. A Novel Spike Distance. Neural Computation, Vol. 13(4), pp. 751-763, 2001.
#' @references S. Schreiber, J.M. Fellous, P.H. Tiesinga, and T.J. Sejnowski. A new correlation-based measure of spike timing reliability, Neurocomputing, Vols. 52-54, pp. 925-931, 2003.
#' @export
#' @examples
#' ##The aftershock data of 26th July 2003 earthquake of M6.2 at the northern Miyagi-Ken Japan.
#' data(Miyagi20030626)
#' ## time longitude latitude depth magnitude 
#' ## split events by 7-hour
#' sMiyagi <- splitMPP(Miyagi20030626,h=60*60*7,scaleMarks=TRUE)$S
#' N <- 10
#' tau <- 0.1
#' sMat <- matrix(0,N,N)
#'   cat("calculating fmetric with tau ",tau,"...")
#'  for(i in 1:(N)){
#'    cat(i," ")
#'    for(j in i:N){
#'      S1 <- sMiyagi[[i]]$time;S2 <- sMiyagi[[j]]$time
#'     sMat[i,j] <- fmetric(S1,S2,tau=tau,M=diag(1,4))
#'    }
#'  }
#'  sMat <- sMat+t(sMat)
#'  tmpd <- diag(sMat) <- diag(sMat)/2
#'  sMat <- sMat/sqrt(outer(tmpd,tmpd))
#' image(sMat)
fmetric <- function(S1,S2,measure="sim",h="laplacian", tau=1,M=NULL,abs.tol=.Machine$double.eps^0.25){
  ## extract information from S1 and S2
  ret <- characterize(S1,S2); T1 <- ret$T1;T2 <- ret$T2;N1 <- ret$N1;N2 <- ret$N2;n.mark <- ret$n.mark
  ## calculate similarity
  if(measure=="sim"){
    ## simple PP
    if(n.mark==0){
      
      if(class(h)=="character"){
        if(h=="laplacian"){
          ##return(as.numeric(0.5*tau*sum(exp(-abs(outer(T1,T2,FUN="-"))/tau))))
          return(as.numeric((1/(4*tau*N1*N2))*sum(exp(-abs(outer(T1,T2,FUN="-"))/tau)))/2)
          ##}
        }else{
          stop("smoother function must be appropriately specified")
        }
      }else if(class(h)=="function"){
        wv <- function(x){
          sum(h(x-T1,tau=tau)*hev(x-T1))*sum(h(x-T2,tau=tau)*hev(x-T2))
        }
        vv <- Vectorize(wv)
        val <- try(integrate(vv,lower=0,upper=Inf,abs.tol=abs.tol)$val,silent=TRUE)
        if(class(val)=="try-error"){
          val <- try(integrate(vv,lower=0,upper=Inf,abs.tol=1e-1,subdivisions=1000)$val,silent=TRUE)
        }
        return(val/(4*N1*N2*tau^2))
      }else{
        stop("smoother function must be appropriately specified")
      }
    }else{  ## MPP
      R1 <- as.matrix(S1[,-1]);R2 <- as.matrix(S2[,-1])
      if(is.null(M)){
        if(n.mark==1){
          tmp <- 1/apply(rbind(R1,R2),2,var);M <- as.matrix(ifelse(is.infinite(tmp),as.matrix(1),(tmp)))
        }else{
          tmp <- 1/apply(rbind(R1,R1),2,var);id <- which(is.infinite(tmp)); if(length(id)){ tmp[id] <- 1}; M <- diag(tmp)
        }
      }
      ## check positive definiteness of precision matrix M
      if(sum(eigen(M)$values <0)){
        stop("precision matrix of smoothing function for marks must be positive definite")
      }
      
      if(class(h)=="character"){
        if(h=="laplacian"){
          if(n.mark==1){
            crossR1R2 <- outer(as.numeric(R1),as.numeric(R2),FUN="-");W <- exp(-(crossR1R2^2)*as.numeric(M/4))
          }else{
            crossR1R2 <- outer(t(R1),t(R2),FUN="-")
            wGauss <- function(x) exp(-mahalanobis(x,center=FALSE,cov=M,inverted=TRUE)/4)
            W <- apply(apply(crossR1R2,c(2,4),diag),c(2,3),wGauss)
          }
          return(as.numeric(sum(W*exp(-abs(outer(T1,T2,FUN="-")/tau)))*sqrt(det(M))/( (2^(n.mark+2))*(pi^(n.mark/2))*tau*N1*N2)))
        }else{
          stop("smoother function must be appropriately specified")
        }
      }else if(class(h)=="function"){
        if(n.mark==1){
          crossR1R2 <- outer(as.numeric(R1),as.numeric(R2),FUN="-");W <- exp(-(crossR1R2^2)*as.numeric(M/4))
        }else{
          crossR1R2 <- outer(t(R1),t(R2),FUN="-")
          wGauss <- function(x) exp(-mahalanobis(x,center=FALSE,cov=M,inverted=TRUE)/4)
          W <- apply(apply(crossR1R2,c(2,4),diag),c(2,3),wGauss)
        }
        wv <- function(x){
          as.numeric((h(x-T1,tau=tau)*hev(x-T1))%*%W%*%(h(x-T2,tau=tau)*hev(x-T2)))
        }
        vv <- Vectorize(wv)
        val <- try(integrate(vv,lower=0,upper=Inf,abs.tol=abs.tol)$val,silent=TRUE)
        if(class(val)=="try-error"){
          val <- try(integrate(vv,lower=0,upper=Inf,abs.tol=1e-1,subdivisions=1000)$val,silent=TRUE)
        }

        val <- val * sqrt(det(M))/( (2^(n.mark+2))*(pi^(n.mark/2))*tau*N1*N2)
        return(val)
      }else{
        stop("smoother function must be appropriately specified")
      }
    }
  }else if(measure=="dist"){  ## calculate distance
    if(n.mark==0){
      if(class(h)=="character"){
        if(h=="laplacian"){
          val <- sqrt(abs(as.numeric(sum(exp(-abs(outer(T1,T1,FUN="-"))/tau)))/(4*tau*N1^2)+as.numeric(sum(exp(-abs(outer(T2,T2,FUN="-"))/tau)))/(4*tau*N2^2)-as.numeric(sum(exp(-abs(outer(T1,T2,FUN="-"))/tau)))/(2*tau*N1*N2)))
          return(val)
        }else{
          stop("smoother function must be appropriately specified")
        }
      }else if(class(h)=="function"){
        wv <- function(x){
          (sum((h(x-T1,tau=tau)*hev(x-T1)))/(2*N1)-sum((h(x-T2,tau=tau)*hev(x-T2)))/(2*N2))^2
        }
        vv <- Vectorize(wv)
        val <- try(integrate(vv,lower=0,upper=Inf,abs.tol=abs.tol)$val,silent=TRUE)
        if(class(val)=="try-error"){
          val <- try(integrate(vv,lower=0,upper=Inf,abs.tol=1e-1,subdivisions=1000)$val,silent=TRUE)
        }
        return(sqrt(val))
      }else{
        stop("smoother function must be appropriately specified")
      }
    }else{  ## Marked PP
      R1 <- as.matrix(S1[,-1]);R2 <- as.matrix(S2[,-1])
      if(is.null(M)){
        if(n.mark==1){
          tmp <- 1/apply(rbind(R1,R2),2,var);M <- as.matrix(ifelse(is.infinite(tmp),as.matrix(1),(tmp)))
        }else{
          tmp <- 1/apply(rbind(R1,R1),2,var);id <- which(is.infinite(tmp)); if(length(id)){ tmp[id] <- 1}; M <- diag(tmp)
        }
      }
      ## check positive definiteness of precision matrix M
      if(sum(eigen(M)$values <0)){
        stop("precision matrix of smoothing function for marks must be positive definite")
      }

      if(n.mark==1){
        crossR1R2 <- outer(as.numeric(R1),as.numeric(R2),FUN="-");W12 <- exp(-(crossR1R2^2)*as.numeric(M/4))
        crossR1 <- outer(as.numeric(R1),as.numeric(R1),FUN="-");W1 <- exp(-(crossR1^2)*as.numeric(M/4))
        crossR2 <- outer(as.numeric(R2),as.numeric(R2),FUN="-");W2 <- exp(-(crossR2^2)*as.numeric(M/4))
      }else{
        crossR1 <- outer(t(R1),t(R1),FUN="-")
        crossR2 <- outer(t(R2),t(R2),FUN="-")
        crossR1R2 <- outer(t(R1),t(R2),FUN="-")
        wGauss <- function(x) exp(-mahalanobis(x,center=FALSE,cov=M,inverted=TRUE)/4)
        W1 <- apply(apply(crossR1,c(2,4),diag),c(2,3),wGauss)
        W2 <- apply(apply(crossR2,c(2,4),diag),c(2,3),wGauss)
        W12 <- apply(apply(crossR1R2,c(2,4),diag),c(2,3),wGauss)
      }
      
      if(class(h)=="character"){
        
        if(h=="laplacian"){
          val <- as.numeric(sum(W1 * exp(-abs(outer(T1,T1,FUN="-"))/(tau))))/(N1^2)+ as.numeric(sum((W2 * exp(-abs(outer(T2,T2,FUN="-"))/(tau)))))/(N2^2) -2*as.numeric(sum((W12 * exp(-abs(outer(T1,T2,FUN="-"))/tau))))/(N1*N2)
          val <- val*sqrt(det(M))/((2^(n.mark+2))*(pi^(n.mark/2))*tau)
          return(sqrt(val))
        }else{
          stop("smoother function must be appropriately specified")
        }
      }else if(class(h)=="function"){
        wv <- function(x){
          (h(x-T1,tau=tau)*hev(x-T1))%*%W1%*%(h(x-T1,tau=tau)*hev(x-T1))/(N1^2)+(h(x-T2,tau=tau)*hev(x-T2))%*%W2%*%(h(x-T2,tau=tau)*hev(x-T2))/(N2^2)-2*(h(x-T1,tau=tau)*hev(x-T1))%*%W12%*%(h(x-T2,tau=tau)*hev(x-T2))/(N1*N2)
        }
        vv <- Vectorize(wv)
        val <- try(integrate(vv,lower=0,upper=Inf,abs.tol=abs.tol)$val,silent=TRUE)
        if(class(val)=="try-error"){
          val <- try(integrate(vv,lower=0,upper=Inf,abs.tol=1e-1,subdivisions=1000)$val,silent=TRUE)
        }
        val <- val * sqrt(det(M))/((2^(n.mark+2))*(pi^(n.mark/2))*tau)
        return(sqrt(val))
      }else{
        stop("smoother function must be appropriately specified")
      }
    }
  }else{
    stop("Measure must be dist or sim")
  }
}
