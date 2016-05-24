#' Metrics for Point Process Realizations Based on Co-occurrence 
#'
#' For comparing two SPP realizations, it is natural to count the number of events which can be considered to be co-occurring. There are two metrics for SPP realizations based on the notion of co-occurrence.
#' The first one proposed by Quian Quiroga et al. (2002) directly counts near-by events. The second counting metric co-occurrence is proposed by Hunter and Milton (2003), which is based on a smoothing function.
#' 
#'
#' \code{coocmetric} computes co-occurrence base metrics for two point process realizations. This function counts the number of events in S1 which is coincided with those in S2, and vice versa.
#'
#'
#' @param S1 marked point process data.
#' @param S2 marked point process data.
#' @param measure \code{"sim"} for similarity and "dist" for distance. Default \code{"sim"}.
#' @param type if \code{"count"}, counting near-by event measure by Quian is computed. If \code{"smooth"}, smoothed counting co-occurrence measure by Hunter and Milton is computed. Default \code{"count"}.
#' @param tau a parameter for filtering function.
#' @param M a precision matrix for filter of marks, i.e., exp( - r' M r) is used for filtering marks. It should be symmetric and positive semi-definite.
#' @return Similarity or distance between two inputs (marked) point process S1 and S2.
#' @author Hideitsu Hino \email{hinohide@@cs.tsukuba.ac.jp}, Ken Takano, Yuki Yoshikawa, and Noboru Murata
#' @references R. Quian Quiroga, T. Kreuz, and P. Grassberger. Event synchronization: a simple and fast method to measure synchronicity and time delay patterns, Physical Review E, Vol. 66(4), 041904, 2002.
#' @references J. D. Hunter and G. Milton. Amplitude and frequency dependence of spike timing: implications for dynamic regulation, Journal of Neurophysiology, Vol. 90, pp. 387-94, 2003.
#' @export
#' @examples
#' ## The aftershock data of 26th July 2003 earthquake of M6.2 at the northern Miyagi-Ken Japan.
#' data(Miyagi20030626)
#' ## time longitude latitude depth magnitude 
#' ## split events by 7-hour
#' sMiyagi <- splitMPP(Miyagi20030626,h=60*60*7,scaleMarks=TRUE)$S
#' N <- 10
#' sMat <- matrix(0,N,N)
#' tau<-0.2
#'   cat("calculating coocmetric(smooth)...")
#'  for(i in 1:(N)){
#'    cat(i," ")
#'    for(j in i:N){
#'      S1 <- sMiyagi[[i]]$time;S2 <- sMiyagi[[j]]$time
#'     sMat[i,j] <- coocmetric(S1,S2,type="smooth",tau=tau,M=diag(1,4))
#'    }
#'  }
#'  sMat <- sMat+t(sMat)
#'  tmpd <- diag(sMat) <- diag(sMat)/2
#'  sMat <- sMat/sqrt(outer(tmpd,tmpd))
#' image(sMat)
coocmetric <- function(S1,S2,measure="sim",type="count",tau=1,M=NULL){
  ## extract information from S1 and S2
  ret <- characterize(S1,S2); T1 <- ret$T1;T2 <- ret$T2;N1 <- ret$N1;N2 <- ret$N2;n.mark <- ret$n.mark
  if(type=="count"){
    ## SPP
    if(n.mark==0){  
      d1 <- d2 <- 0
      grid <- expand.grid(1:N1,1:N2)
      apply(grid,1,function(x){
        k1 <- x[1]
        k2 <- x[2]
        g <- c(T1[k1+1]-T1[k1],T1[k1]-T1[k1-1],T2[k2+1]-T2[k2],T2[k2]-T2[k2-1])
        g <- min(na.exclude(g))/2
        
        id1 <- T1[k1]-T2[k2]
        id2 <- T2[k2]-T1[k1]
        
        if( (0 < id1) & (id1 < g)){
          d1 <<- d1 + 1
        }else if(id1==0){
          d1 <<- d1 + 1/2
        }else{
        }
        
        if((0 < id2) & (id2 < g)){
          d2 <<- d2 + 1
        }else if(id2==0){
          d2 <<- d2 + 1/2
        }else{
        }
        
      })
      val <- (d1 + d2) / sqrt(N1*N2)
      if(measure=="sim"){  ## calculate similarity
        return(val)
      }else if(measure=="dist"){  ## calculate distance
        return(1-val)
      }else{
        stop("Measure must be dist or sim")
      }
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
         crossR1R2 <- outer(as.numeric(R1),as.numeric(R2),FUN="-");W <- exp(-(crossR1R2^2)*as.numeric(M))
       }else{
         crossR1R2 <- outer(t(R1),t(R2),FUN="-")
         wGauss <- function(x) exp(-mahalanobis(x,center=FALSE,cov=M,inverted=TRUE)/2)
         W <- apply(apply(crossR1R2,c(2,4),diag),c(2,3),wGauss)
       }
      
      d1 <- d2 <- 0
      grid <- expand.grid(1:N1,1:N2)
      apply(grid,1,function(x){
        k1 <- x[1]
        k2 <- x[2]
        g <- c(T1[k1+1]-T1[k1],T1[k1]-T1[k1-1],T2[k2+1]-T2[k2],T2[k2]-T2[k2-1])
        g <- min(na.exclude(g))/2
        
        id1 <- T1[k1]-T2[k2]
        id2 <- T2[k2]-T1[k1]
        
        if( (0 < id1) & (id1 < g)){
          d1 <<- d1 + W[k1,k2]
        }else if(id1==0){
          d1 <<- d1 + 0.5*W[k1,k2]
        }else{
        }
        
        if((0 < id2) & (id2 < g)){
          d2 <<- d2 + W[k1,k2]
        }else if(id2==0){
          d2 <<- d2 + 0.5*W[k1,k2]
        }else{
        }
      })
      val <- (d1 + d2) / sqrt(N1*N2)
      if(measure=="sim"){  ## calculate similarity
        return(val)
      }else if(measure=="dist"){  ## calculate distance
        return(1-val)
      }else{
        stop("Measure must be dist or sim")
      }
    }
  }else if(type=="smooth"){
    x <- as.numeric(c(T1,T2))
    names(x)<- c(rep(1,N1),rep(2,N2))
    x1ind <- which(names(x)==1)
    x2ind <- which(names(x)==2)
    names(x) <- NULL
    xid <- c(x1ind,seq(1,N2))
    
    dx1 <- dx2 <- 0
    
    if(n.mark!=0){
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
    }
      
    for(i in 1:length(x)){
      if(i %in% x1ind){
        nnID <- which.min(abs(x[x2ind]-x[i]))
        dx1 <- dx1 + exp(-min(abs(x[x2ind]-x[i]))/tau)*ifelse(n.mark==0,1, exp(-mahalanobis( (R1[xid[i],]-R2[xid[nnID],]),center=FALSE,cov=M,inverted=TRUE)))
      }else{
        nnID <- which.min(abs(x[x1ind]-x[i]))
        dx2 <- dx2 + exp(-min(abs(x[x1ind]-x[i]))/tau)*ifelse(n.mark==0,1, exp(-mahalanobis( (R2[xid[i],]-R1[xid[nnID],]),center=FALSE,cov=M,inverted=TRUE)))
      }
    }
    val <- 0.5*(dx1/N1 + dx2/N2)
    if(measure=="sim"){  ## calculate similarity
      return(val)
    }else if(measure=="dist"){  ## calculate distance
      return(1-val)
    }else{
      stop("Measure must be dist or sim")
    }
  }else{
    stop("Type must be count or smooth")
  }
}
