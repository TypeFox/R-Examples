#' Compute Intensity Inner Product Metrics
#'
#'
#' For the analysis of point process, intensity function plays a central roll. Paiva et al. (2009) proposed to use the intensity function for defining the inner product between point process realizations. 
#'
#'
#' \code{iipmetric} computes intensity inner product metric. Intensity function for the point process realization is estimated by kernel density estimator. This function adopts Gaussian kernels for the sake of computational efficiency.
#' @param S1 marked point process data.
#' @param S2 marked point process data.
#' @param measure \code{"sim"} for similarity and \code{"dist"} for distance. Default \code{"sim"}.
#' @param tau a parameter for filtering function.
#' @param M a precision matrix for filter of marks, i.e., exp( - r' M r) is used for filtering marks. It should be symmetric and positive semi-definite.
#' @return Similarity or distance between two inputs (marked) point process S1 and S2.
#' @author Hideitsu Hino \email{hinohide@@cs.tsukuba.ac.jp}, Ken Takano, Yuki Yoshikawa, and Noboru Murata
#' @references A.R.C. Paiva, I. Park, and J.C. Principe. A reproducing kernel Hilbert space framework for spike train signal processing, Neural Computation, Vol. 21(2), pp. 424-449, 2009.
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
#'   cat("calculating intensity inner product...")
#'  for(i in 1:(N)){
#'    cat(i," ")
#'    for(j in i:N){
#'      S1 <- sMiyagi[[i]]$time;S2 <- sMiyagi[[j]]$time
#'     sMat[i,j] <- iipmetric(S1,S2,tau=tau,M=diag(1,4))
#'    }
#'  }
#'  sMat <- sMat+t(sMat)
#'  tmpd <- diag(sMat) <- diag(sMat)/2
#'  sMat <- sMat/sqrt(outer(tmpd,tmpd))
#' image(sMat)
iipmetric <- function(S1,S2,measure="sim",tau=1,M=NULL){
  ## extract information from S1 and S2
  ret <- characterize(S1,S2); T1 <- ret$T1;T2 <- ret$T2;N1 <- ret$N1;N2 <- ret$N2;n.mark <- ret$n.mark
  ## calculate similarity
  if(measure=="sim"){
    ## ground intensity function
    gi <- as.numeric((exp(-(outer(T1,T2,FUN="-")^2)/(4*tau^2))))/(4*sqrt(pi)*tau*N1*N2)
    if(n.mark==0){
      return(sum(gi))
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
      if(n.mark==1){
        crossR1R2 <- outer(as.numeric(R1),as.numeric(R2),FUN="-");W <- exp(-(crossR1R2^2)*as.numeric(M/4))
      }else{
        crossR1R2 <- outer(t(R1),t(R2),FUN="-")
        wGauss <- function(x) exp(-mahalanobis(x,center=FALSE,cov=M,inverted=TRUE)/4)
        ## density for marks
        W <- apply(apply(crossR1R2,c(2,4),diag),c(2,3),wGauss)
      }
      return(sum(gi*W)/((pi^(n.mark/2))))
    } 
  }else if(measure=="dist"){ ## calculate distance
    gij <- as.numeric((exp(-(outer(T1,T2,FUN="-")^2)/(4*tau))))/(4*sqrt(pi)*tau*N1*N2)
    gii <- as.numeric((exp(-(outer(T1,T1,FUN="-")^2)/(4*tau))))/(4*sqrt(pi)*tau*N1^2)
    gjj <- as.numeric((exp(-(outer(T2,T2,FUN="-")^2)/(4*tau))))/(4*sqrt(pi)*tau*N2^2)
    if(n.mark==0){
      return(sqrt(sum(gii)+sum(gjj)-2*sum(gij)))
    }else{  ## MPP
      R1 <- as.matrix(S1[,-1]);R2 <- as.matrix(S2[,-1])
      if(is.null(M)){
        if(n.mark==1){
          tmp <- 1/apply(rbind(R1,R2),2,var);M <- as.matrix(ifelse(is.infinite(tmp),as.matrix(1),(tmp)))
        }else{
          tmp <- 1/apply(rbind(R1,R1),2,var);id <- which(is.infinite(tmp)); if(length(id)){ tmp[id] <- 1}; M <- diag(tmp)
          ##M <- diag(1/apply(rbind(R1,R2),2,var))
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
        crossR1R2 <- outer(t(R1),t(R2),FUN="-")
        crossR1 <- outer(t(R1),t(R1),FUN="-")
        crossR2 <- outer(t(R2),t(R2),FUN="-")
        wGauss <- function(x) exp(-mahalanobis(x,center=FALSE,cov=M,inverted=TRUE)/4)
        W1 <- apply(apply(crossR1,c(2,4),diag),c(2,3),wGauss)
        W2 <- apply(apply(crossR2,c(2,4),diag),c(2,3),wGauss)
        W12 <- apply(apply(crossR1R2,c(2,4),diag),c(2,3),wGauss)
      }
      val <- (as.numeric(sum(W1 * exp(-(outer(T1,T1,FUN="-")^2)/(2*tau)))/(N1^2))+
              as.numeric(sum((W2 * exp(-(outer(T2,T2,FUN="-")^2)/(2*tau))))/(N2^2))+
              -2*as.numeric(sum((W12 * exp(-(outer(T1,T2,FUN="-")^2)/(2*tau)))))/(N1*N2))
      val <- val/((pi^((n.mark+1)/2))*tau*sqrt(det(M)))
      return(sqrt(val))
    }
  }else{
    stop("Measure must be dist or sim")
  }
}
