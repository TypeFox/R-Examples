##### I am using a method similar to cartsplit for survival since I am using the same L2 loss 
# function just with different weights. For IPCW 
## and Brier when there is only 1 cutpoint, we can directly use the cartsplit function. 
# However, when have multiple cutpoints, I am attempting
## to use the same code as in cartsplit, but I'm averaging goodness over all the cutpoints 
# to choose the best split point. This concept was developed
## in the interest of keeping the running time reasonable. The approach I took over the summer 
# results in a very high running time, so I am hoping
## that this will be an improvement.
survival.split <- function(psi, y, wt,  x.split, minsplit, minbuck, real.num, opts,is.num) {


## Simple case where the weights are constant
if(opts$loss.fx=="IPCW") {
   
                
 return( cartsplit(psi=psi, y=y[,1], wt=wt, x.split=x.split,minsplit=minsplit,  minbuck=minbuck, real.num=real.num, opts=opts,is.num=is.num) )


}else if(opts$loss.fx=="Brier" && length(opts$brier.vec)==1){
   
	y=as.numeric(y[,1]>opts$brier.vec[1])
    
    return( cartsplit(psi=psi, y=y, wt=wt, x.split=x.split, minsplit=minsplit, minbuck=minbuck, real.num=real.num, opts=opts,is.num=is.num) )
## more complicated case where we have changing y and wt values depending on which cutpoint we use
}else if(opts$loss.fx=="Brier" && length(opts$brier.vec)>1){  
     
    wt.Brier=wt[,1]
    
     y.Brier=as.numeric(y[,1]>opts$brier.vec[1])

     
     if(sum(is.na(x.split))>0){missing="yes"}else{missing="no"}
  
  
  
  
  
  ## psi tells which observations are in node that we are trying to split
  
  
 
 
  STOP <- 0
  x.s <- x.split[order(x.split)]
  wt.s <- wt.Brier[order(x.split)]
  y.s <- y.Brier[order(x.split)]
  psi.s <- psi[order(x.split)]
  x.s.k <- x.s[psi.s == 1]
  wt.s.k <- wt.s[psi.s == 1]
  y.s.k <- y.s[psi.s == 1]

  
  
  
  
                                        # if there are missing values, impute the x to be x.impute using the mean or
                                        #  mode. then update x.s.k, y.s.k, wt.s.k to be the non-missing values since
                                        #  that is what will be used to create the splits.
  if(missing=="yes"){
    
    x.original <- x.s.k
    x.impute<-x.split
    
    x.impute[which(is.na(x.split))] <- ifelse(is.num==1,mean(x.s.k, na.rm=TRUE),as.numeric(names(which.max(table(x.s.k)))))
    
    toSave=which(!is.na(x.s.k))
    y.s.k <- y.s.k[toSave]
    wt.s.k <- wt.s.k[toSave]
    x.s.k <- x.s.k[toSave]
    
    
    
  }
  
  
  
  
  
  y.s.k <- y.s.k- sum(y.s.k * wt.s.k) / sum(wt.s.k)
  cant.split <- FALSE
  if(length(unique(x.s.k)) == 1 | length(y.s.k)<minsplit | length(which(wt.s.k!=0))<=1) {
    
    STOP <- 1
  } else if(abs(minsplit - max(table(x.s.k))) < minbuck &&
            ((length(x.s.k) - max(table(x.s.k))) < minbuck)) {
    STOP <- 1
  } else {  
    
	if (length(unique(y.s.k))!=1){
    wtsum <- tapply(wt.s.k, x.s.k, sum)
    wtsum<-wtsum[!is.na(wtsum)]
    ysum  <- tapply(y.s.k * wt.s.k, x.s.k, sum)
    ysum<-ysum[!is.na(ysum)]
    means <- ysum / wtsum

    
    n <- length(ysum)
    temp <- cumsum(ysum)[-n]
    left.wt  <- cumsum(wtsum)[-n]
    right.wt <- sum(wt.s.k) - left.wt
    lmean <- temp/left.wt
    rmean <- -temp/right.wt
    
    goodness <- opts$IBS.wt[1] * (left.wt * lmean^2 + right.wt * rmean^2)  / sum(wt.s.k * y.s.k^2)
  }else{goodness <- NULL}
      
    ##### Need to take into account other brier.vec cutpoints so we loop over all of them and update the weights and y values
    for (k in 2:ncol(wt)){
    
      ## recompute wt's and y's and get wt.s.k and y.s.k
  
      wt.Brier=wt[,k]
      y.Brier=as.numeric(y[,1]>opts$brier.vec[k])
     
      wt.s <- wt.Brier[order(x.split)]
      y.s <- y.Brier[order(x.split)]
      wt.s.k <- wt.s[psi.s == 1]
      y.s.k <- y.s[psi.s == 1]
  
      if(missing=="yes"){
       y.s.k <- y.s.k[toSave]
       wt.s.k <- wt.s.k[toSave]
       }
       
       y.s.k <- y.s.k- sum(y.s.k * wt.s.k) / sum(wt.s.k)
       ## we only want to include this in goodness if there are multiple y.s.k values
       if(length(unique(y.s.k))!=1){
  
           ## calculate a new goodness called next.goodness using same method
           wtsum <- tapply(wt.s.k, x.s.k, sum)
           wtsum<-wtsum[!is.na(wtsum)]
           ysum  <- tapply(y.s.k * wt.s.k, x.s.k, sum)
           ysum<-ysum[!is.na(ysum)]
           means <- ysum / wtsum
    
    
           n <- length(ysum)
           temp <- cumsum(ysum)[-n]
           left.wt  <- cumsum(wtsum)[-n]
           right.wt <- sum(wt.s.k) - left.wt
           lmean <- temp/left.wt
    
           rmean <- -temp/right.wt
           next.goodness <-(left.wt * lmean^2 + right.wt * rmean^2) / sum(wt.s.k * y.s.k^2)
           ## combine all goodnesses into a matrix
           goodness=rbind(goodness,opts$IBS.wt[k] * next.goodness)
        }
  
    } 
    
    ## take a mean over the columns of this matrix
    if(!is.vector(goodness) & !is.null(goodness)){
	goodness=apply(goodness,2,mean,na.rm=TRUE)
	}
	 ### added 10/23/11 - this checks if the goodness values are all NA and if so, it tells us that we can't do a reasonable split on this node so we set STOP=1.
     ### I'm pretty sure that this scenario corresponds to not having any group 2 observations in the node for any of the time points. And if we only have group 1 and 3 observations,
      ### then in a sense isn't the node already pure?

    if(is.null(goodness)) {
         STOP<-1
      }else if (sum(is.na(goodness))==length(goodness)){
         STOP <-1
      }else{


    
    best.ind <- central.split(which.max(goodness),n)
    
    max.g <- max(goodness,na.rm=T)
    in.set <- as.numeric(names(goodness)[best.ind])
    
    
    if(missing=="no") {in.x <- (x.split <= (in.set + 1e-10))}
    if(missing=="yes") {in.x <- (x.impute <= (in.set + 1e-10))}
    num.in.set <- sum(in.x & psi==1 & wt.Brier>0)
    num.not.in.set <- length(which(wt.s.k>0)) - num.in.set
    
    
    if((num.in.set >= minbuck) &&
       (num.not.in.set >= minbuck)) {
      
      if(missing=="no") {in.x <- (x.split <= (in.set + 1e-10))}
      if(missing=="yes") {in.x <- (x.impute <= (in.set + 1e-10))}
      
      new.lt <- as.numeric(in.x)
      new.gt <- 1 - new.lt
      new.lt <- ifelse(psi == 0, 0, new.lt)
      new.gt <- ifelse(psi == 0, 0, new.gt)
       val <- in.set + 1e-10
      cant.split <- FALSE
    } else {
      ## so we don't split a variable resulting in
      ## a vector with less than minbuck observations
      n.chg <- 2
      
      ## sequentially try goodness values until we have a split
      ## with enough observations in each basis functions
      ## added in is.na(max.g) because it is possible that some goodness values are NA if first or last weights are 0.
      while(is.na(max.g) | !(num.in.set >= minbuck &&
              num.not.in.set >= minbuck)) {
        
        if(n.chg >= n) {
          STOP <- 1
          break
        }
        max.g <- goodness[order(goodness)][n - n.chg]
        n.chg <- n.chg + 1
        best.ind <- central.split(which(goodness == max.g),n)
        in.set <- as.numeric( names(goodness)[best.ind])
        if(missing=="no"){ in.x <- (x.split <= (in.set+1e-10) )}
        if(missing=="yes"){in.x <- (x.impute <= (in.set + 1e-10))}
        num.in.set <- sum(in.x & psi==1 & wt.Brier>0)
        num.not.in.set <- length(which(wt.s.k>0)) - num.in.set
        
        
        
      } # end while statement
      
      if(missing=="no"){in.x <- (x.split <= (in.set + 1e-10))}
      if(missing=="yes"){in.x <- (x.impute <= (in.set + 1e-10))}
      new.lt <- as.numeric(in.x)
      new.gt <- 1 - new.lt
      new.lt <- ifelse(psi == 0, 0, new.lt)
      new.gt <- ifelse(psi == 0, 0, new.gt)
      val <-in.set + 1e-10
      cant.split <- FALSE
      
    } # end else
    } # end else
	
  } # end main code
  if(STOP )  {
    new.lt <- psi
    new.gt <- rep(0, length(psi))
    cant.split <- TRUE
    best.ind <- 1
    val <- NA
    max.g <- NA
  }
  
  return(list(new.lt=new.lt, new.gt=new.gt,
              val=val, cant.split=cant.split,
              max.g=max.g))
}  # ends the multiple brier vector scenario

    
   
  
} #end outer function



