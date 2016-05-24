#' Compute Inter Event Interval-based Metric Between Marked Point Processes
#'
#' 
#' This metric considers inter event interval for point processes. 
#'
#'
#' \code{iei} computes inter event interval-based measure between MPP realizations. iei for simple point process does not have any tuning parameter, which can be a desirable property for data analysis. However, it's computational cost is relatively higher than other metrics.
#' 
#' @param S1 marked point process data.
#' @param S2 marked point process data.
#' @param measure \code{"sim"} for similarity and \code{"dist"} for distance. Default \code{"sim"}.
#' @param M a precision matrix for filter of marks, i.e., exp( - r' M r) is used for filtering marks. It should be symmetric and positive semi-definite.
#' @param window.length  width of the window used for splitting the original MPP.\cr
#' If not provided, \code{max(max(S1$time,S2$time) - min(S1$time,S2$time))} is used.
#' @param variant  choose from two variants "spike-weighted" or "time-weighted".\cr Default \code{"spike"}, which is computationally efficient than \code{"time"}. See the reference for details.
#' @param abs.tol absolute tolerance for numerical integration.
#' @return Similarity or distance between two inputs (marked) point process S1 and S2.
#' @author Hideitsu Hino \email{hinohide@@cs.tsukuba.ac.jp}, Ken Takano, Yuki Yoshikawa, and Noboru Murata
#' @export
#' @references T. Kreuz, J.S. Haas, A. Morelli, H.D.I. Abarbanel, and A. Politi. Measuring spike train synchrony, Journal of Neuroscience Methods, Vol. 165(1), pp. 151-161, 2007.
#' @examples
#' ##The aftershock data of 26th July 2003 earthquake of M6.2 at the northern Miyagi-Ken Japan.
#' data(Miyagi20030626)
#' ## time longitude latitude depth magnitude 
#' ## split events by 7-hour
#' sMiyagi <- splitMPP(Miyagi20030626,h=60*60*7,scaleMarks=TRUE)$S
#' N <- 5
#' sMat <- matrix(0,N,N)
#'   cat("calculating intensity inner product...")
#'  for(i in 1:(N)){
#'    cat(i," ")
#'    for(j in i:N){
#'      S1 <- sMiyagi[[i]]$time;S2 <- sMiyagi[[j]]$time
#'     sMat[i,j] <- ieimetric(S1,S2,M=diag(1,4))
#'    }
#'  }
#'  sMat <- sMat+t(sMat)
#'  tmpd <- diag(sMat) <- diag(sMat)/2
#'  sMat <- sMat/sqrt(outer(tmpd,tmpd))
#' image(sMat)
ieimetric <- function(S1,S2,measure="sim",M=NULL,window.length=NULL,variant="spike",abs.tol=.Machine$double.eps^0.25){
  ## extract information from S1 and S2
  ret <- characterize(S1,S2); T1 <- ret$T1;T2 <- ret$T2;N1 <- ret$N1;N2 <- ret$N2;n.mark <- ret$n.mark
  if(setequal(T1,T2) && n.mark==0){return(ifelse(measure=="sim",1,0))}  
  ## preprocessing
  if(is.null(window.length)){
    window.length <- max(max(max(c(T1,T2)) - min(c(T1,T2))),1)
  }
  
  ##--------------------------------------------------------------
  ## iei: returns the iei value at time t for pp data S
  ##--------------------------------------------------------------
  v <- function(t,S){
    if(t <= S[1]){
      ieivalue <- S[1]
    }else if(t>=S[length(S)]){
      ieivalue <- window.length - S[length(S)]
    }else{
      tmp1 <- subset(S, S<t)
      maxid <- which.max(tmp1)
      tmp2 <- subset(S, S>t)
      minid <- which.min(tmp2)
      ieivalue <- tmp2[minid]-tmp1[maxid]
    }
    return(ieivalue)
  }

  ## intermediate function
  imd <- function(x,y,t){
    xiei <- v(S=x,t=t)
    yiei <- v(S=y,t=t)
    if(xiei==0 | yiei==0){return(1)}  ## 2014/11/05 fixed by Ken Takano
    val <- min(xiei,yiei)/max(xiei,yiei)
    return(val)
  }
  if(n.mark==0){  ## SPP
    if(measure=="sim"){
      if(variant=="time"){
        wimd <- function(t){      imd(T1,T2,t)    }
        vimd <- Vectorize(wimd)
        val <- try(integrate(vimd,lower=0,upper=window.length,abs.tol=abs.tol)$val,silent=TRUE)
        if(class(val)=="try-error"){
          val <- try(integrate(vimd,lower=0,upper=window.length,abs.tol=1e-1,subdivisions=1000)$val,silent=TRUE)
        }
        }else if(variant=="spike"){
          eventlist <- unique(c(T1,T2))
          val <- 0
          for(t in eventlist){
            val <- val+imd(T1,T2,t)
          }
        }else{
          stop("var must be either spike or time")
        }
        return(val/window.length)
    }else if(measure=="dist"){
      if(variant=="time"){
        wimd <- function(t){      abs(1-imd(T1,T2,t))    }
        vimd <- Vectorize(wimd)
        val <- try(integrate(vimd,lower=0,upper=window.length,abs.tol=abs.tol)$val,silent=TRUE)
        if(class(val)=="try-error"){
          val <- try(integrate(vimd,lower=0,upper=window.length,abs.tol=1e-1,subdivisions=1000)$val,silent=TRUE)
        }
      }else if(variant=="spike"){
        eventlist <- unique(c(T1,T2))
        val <- 0
        for(t in eventlist){
          val <- val+(1-imd(T1,T2,t))
        }
      }else{
        stop("var must be either spike or time")
      }
      return(val/window.length)
    }else{
      stop("Measure must be dist or sim")
    }
  }else{  ## MPP
    R1 <- as.matrix(S1[,-1]);R2 <- as.matrix(S2[,-1])
    ## check positive definiteness of precision matrix M
    if(is.null(M)){
      if(n.mark==1){
        tmp <- 1/apply(rbind(R1,R2),2,var);M <- as.matrix(ifelse(is.infinite(tmp),as.matrix(1),(tmp)))
        
      }else{
        tmp <- 1/apply(rbind(R1,R1),2,var);id <- which(is.infinite(tmp)); if(length(id)){ tmp[id] <- 1}; M <- diag(tmp)
        ##M <- diag(1/apply(rbind(R1,R2),2,var))
      }
    }
    if(sum(eigen(M)$values <0)){
      stop("precision matrix of smoothing function for marks must be positive definite")
    }
    qul <- function(t){
      if(t <= T1[1]){
        mw1 <- exp(-mahalanobis(R1[1,],center=FALSE,cov=M,inverted=TRUE))
      }else if(t>=T1[length(T1)]){
        mw1 <- exp(-mahalanobis(R1[N1,],center=FALSE,cov=M,inverted=TRUE))
      }else{
        tmp1 <- subset(T1, T1<t)
        maxid <- which.max(tmp1)
        tmp2 <- subset(T1, T1>t)
        minid <- which.min(tmp2)
        mw1 <- exp(-mahalanobis(R1[maxid,],center=FALSE,cov=M,inverted=TRUE))+exp(-mahalanobis(R1[minid,],center=FALSE,cov=M,inverted=TRUE))
      }
      if(t <= T2[1]){
        mw2 <- exp(-mahalanobis(R2[1,],center=FALSE,cov=M,inverted=TRUE))
      }else if(t>=T2[length(T2)]){
        mw2 <- exp(-mahalanobis(R2[N2,],center=FALSE,cov=M,inverted=TRUE))
      }else{
        tmp1 <- subset(T2, T2<t)
        maxid <- which.max(tmp1)
        tmp2 <- subset(T2, T2>t)
        minid <- which.min(tmp2)
        mw2 <- exp(-mahalanobis(R2[maxid,],center=FALSE,cov=M,inverted=TRUE))+exp(-mahalanobis(R2[minid,],center=FALSE,cov=M,inverted=TRUE))
      }
      return(0.5*(mw1+mw2))
    }
    if(measure=="sim"){
      if(variant=="time"){
        wimd <- function(t){ imd(T1,T2,t)*qul(t) }
        vimd <- Vectorize(wimd)
        val <- try(integrate(vimd,lower=0,upper=window.length,abs.tol=abs.tol)$val,silent=TRUE)
        if(class(val)=="try-error"){
          val <- try(integrate(vimd,lower=0,upper=window.length,abs.tol=1e-1,subdivisions=1000)$val,silent=TRUE)
        }
      }else if(variant=="spike"){
        eventlist <- unique(c(T1,T2))
        val <- 0
        for(t in eventlist){
          val <- val+imd(T1,T2,t)*qul(t)
        }
      }else{
        stop("var must be either spike or time")
      }
      return(val/window.length)
    }else if(measure=="dist"){
      if(variant=="time"){
        wimd <- function(t){      abs(1-imd(T1,T2,t)*qul(t))    }
        vimd <- Vectorize(wimd)
        val <- try(integrate(vimd,lower=0,upper=window.length,abs.tol=abs.tol)$val,silent=TRUE)
        if(class(val)=="try-error"){
          val <- try(integrate(vimd,lower=0,upper=window.length,abs.tol=1e-1,subdivisions=1000)$val,silent=TRUE)
        }
      }else if(variant=="spike"){
        eventlist <- unique(c(T1,T2))
        val <- 0
        for(t in eventlist){
          val <- val+(1-imd(T1,T2,t)*qul(t))
        }
      }else{
        stop("var must be either spike or time")
       }
      return(val/window.length)
    }else{
      stop("Measure must be dist or sim")
    }
  }
}
