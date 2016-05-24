#' Split MPP Data by Sliding Time Window
#'
#' This function splits a point process realization into a list of splitted point process realizations with length \code{h}.
#'
#' \code{splitMPP} splits a point process realization into a list of splitted point process realizations with length \code{h}.
#'
#'
#' @param mppdata marked point process in data.frame composed of "time, mark1, mark2, ...".
#' @param h width of the time window. Default is set to h=60*60*48, which is two days when $time is recorded in second. This is suitable for a special seismic data only.
#' @param ol length of overlap for the sliding window. Default 0.
#' @param TimeOrigin logical. If \code{TRUE}, the beginning of the window is assumed to be the origin of time. Default \code{TRUE}.
#' @param scaleMarks logical. If \code{TRUE}, marks (except time) are normalized to have unit variance in whole time series (not in individual windows). Default \code{FALSE}.
#' @param scaleWindow logical. If \code{TRUE}, time interval (window.length) is normalized to one.
#' @param MarkCenter vector for specifying the center of the mark. Use when there are relative center point such as the main shock of the earthquake. Default \code{NULL}.
#' @param MarkCenterID vector for specifying the elements of center of the mark. 
#' @export
#' @examples
#' ##The aftershock data of 26th July 2003 earthquake of M6.2 at the northern Miyagi-Ken Japan.
#' data(Miyagi20030626)
#' ## time longitude latitude depth magnitude 
#' ## split events by 5-hours
#' sMiyagi <- splitMPP(Miyagi20030626,h=60*60*5,scaleMarks=TRUE)
splitMPP <- function(mppdata,h=60*60*48,ol=NULL,TimeOrigin=TRUE,scaleMarks=FALSE,scaleWindow=TRUE,MarkCenter=NULL,MarkCenterID=NULL){
  if(scaleMarks){
    X <- scale(mppdata[,-1])
    mppMeans <- attr(X,"scaled:center")
    mppSDs <- attr(X,"scaled:scale")
    mppdata[,-1] <- X;rm(X)
  }else if(!is.null(MarkCenter)){
    X <- mppdata[,-1]
    mppSDs <- NULL
    if(is.null(MarkCenterID)){
      mppMeans <- MarkCenter
      X <- as.matrix(X)-t(matrix(rep(as.numeric(MarkCenter),time=dim(X)[1]),,dim(X)[1]))
      mppdata[,-1] <- data.frame(X)
    }else{
      for(j in MarkCenterID){
        X[,j] <- X[,j]-as.numeric(MarkCenter[j])
      }
      mppdata[,-1] <- data.frame(X)
      mppMeans <- numeric(dim(mppdata)[2]-1)
      mppMeans[MarkCenterID] <- MarkCenter[MarkCenterID]
    }
  }else{
    mppMeans <- NULL;mppSDs <- NULL
  }
  time.seq <- mppdata$time

  ## the dimension of marks
  markdim <- dim(mppdata)[2]-1

  ## mark names
  marknames <- paste("mark",seq(1,markdim),sep="")
  S <- list()
  noEventID <- NULL
  if( is.null(ol) || ol==0  ){
    ol <- h
  }
  tics.start <- seq(time.seq[1],time.seq[(length(time.seq)-1)],by=ol)
  tics.end <- seq(time.seq[1]+h,time.seq[(length(time.seq)-1)]+h,by=ol)
  for(i in 1:(length(tics.start))){
    id <- which(((time.seq > tics.start[i])*(time.seq < tics.end[i]))==1)
    S[[length(S)+1]] <- mppdata[id,]
    
    ## when there is no event in the window, set the event occuring time to "-1" and marks to "0"
    ## treatment of these windows is left for the user
    if(length(S[[length(S)]][,1])==0){
      noEventID <- c(noEventID,i)
      S[[length(S)]] <- c(-1,rep(0,markdim))
      names(S[[length(S)]]) <- c("time",marknames); S[[length(S)]] <- data.frame(t(S[[length(S)]]))
    }else if(TimeOrigin){   ## if TimeOrigin is TRUE, the event time is recorded from the beggining of the window
      S[[length(S)]]$time <- S[[length(S)]]$time-tics.start[i]
    }
    
    if(scaleWindow){
      S[[length(S)]]$time <- S[[length(S)]]$time/h
    }
  }
  
  names(S) <- tics.start
  return(list(S=S, mppMeans=mppMeans, mppSDs=mppSDs,noEventID=noEventID))
}
