#' Binarize the drug target profile data
#' 
#' A function for binarizing the drug target profile data.
#' 
#' @param profile a matrix with non-binary entries. The rows are drugs and the columns are targets.
#' @param method a string to specify the methods used for binarizing the data. When it is "universal", an
#' universal threshold is used. In such case, another parameter threshold can only be one of "100nM", "1000nM",
#' and "10000nM". When it is "drug-specific", the threshold used for binarization depends on each drug and
#' the parameter threshold can be only one of "10fold", "50fold", and "100fold".
#' @param threshold a string to specify the threshold.
#' @return A matrix contains the binarized drug target data. 
#' 
#' @author Jing Tang \email{jing.tang@@helsinki.fi} 
#' @examples 
#' data(davis)
#' profile<-binarizeDrugTargets(davis, method="drug-specific", threshold="50fold")
#'
binarizeDrugTargets<-function(profile, method="universal", threshold="100nM"){
  # find the NA entries
  nas <- which(is.na(profile)==TRUE, arr.ind=TRUE)
  if(method=="universal"){
    # find entries over threshold -> 0
    if(threshold=="100nM"){
      zeros <- which(profile>=100, arr.ind=TRUE)
      ones <- which(profile<100, arr.ind=TRUE)
      profile[zeros] <- 0
      profile[ones] <- 1
    }else if(threshold=="1000nM"){
      zeros <- which(profile>=1000, arr.ind=TRUE)
      ones <- which(profile<1000, arr.ind=TRUE)
      profile[zeros] <- 0
      profile[ones] <- 1
    }else if(threshold=="10000nM"){
      zeros <- which(profile>=10000, arr.ind=TRUE)
      ones <- which(profile<10000, arr.ind=TRUE)
      profile[zeros] <- 0
      profile[ones] <- 1
    }else{stop("undefined threshold for the universal method.")}
  }else if(method=="drug-specific"){
    if(threshold=="10fold"){
      for(i in seq_len(nrow(profile))){
        # find the min kd
        min_kd <- min(profile[i,], na.rm=TRUE)
        # threshold
        cut_off <- min_kd*10
        if(cut_off>10000) cut_off <- 10000
        ones <- which(profile[i,]<=cut_off, arr.ind=TRUE)
        zeros <- which(profile[i,]>cut_off, arr.ind=TRUE)
        profile[i,ones] <- 1
        profile[i, zeros] <- 0
      }
    }else if(threshold=="50fold"){
      for(i in seq_len(nrow(profile))){
        # find the min kd
        min_kd <- min(profile[i,], na.rm=TRUE)
        # threshold
        cut_off <- min_kd*50
        if(cut_off>10000) cut_off <- 10000
        ones <- which(profile[i,]<=cut_off, arr.ind=TRUE)
        zeros <- which(profile[i,]>cut_off, arr.ind=TRUE)
        profile[i,ones] <- 1
        profile[i, zeros] <- 0
      }
    }else if(threshold=="100fold"){
      for(i in seq_len(nrow(profile))){
        # find the min kd
        min_kd <- min(profile[i,], na.rm=TRUE)
        # threshold
        cut_off <- min_kd*100
        if(cut_off>10000) cut_off <- 10000
        ones <- which(profile[i,]<=cut_off, arr.ind=TRUE)
        zeros <- which(profile[i,]>cut_off, arr.ind=TRUE)
        profile[i,ones] <- 1
        profile[i, zeros] <- 0
      }
    }else {stop("undefined threshold for the drug-specific method.")}
  }else {stop("undefined method.")}
  profile[nas] <- 0
  return(profile)
}