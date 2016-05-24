#########################################################################
###                                                                     #
###  Cartsplit is a function which takes a basis function (indicator    #
###        variable) and returns the best split of that basis function  #
###        into two disjoint sets.                                      #
###  The function is based on code from rpart's anova splitting function#
###                                                                     #
###  Cartsplit requires the following parameters:                       #
###    Psi = indicator function                                         #
###    Psi(Y) = outcome variable of those observations included in Psi  #
###    Psi(X) = variable Psi was built on of those obs included in Psi  #
###    Psi(Wt) = weights of those obs included in Psi                   #
###                                                                     #
###  Cartsplit returns a list of two indicator variables I(Psi(X)<c)    #
###    and I(Psi(X)>c)                                                  #
###                                                                     #
###  Cartsplit is called by addon() and subst()                         #
#########################################################################


##############################################################################
###  Presently we don't allow a cut to occur if it will result in one vector #
###  of one and the rest in the other vector. Instead we pick the next best  #
###  split. However, this is hardcoded - so user can not change - might want #
###  to change in the future                                                 #
##############################################################################

split.fx <- function(psi, y, wt,  x.split, minsplit, minbuck,  real.num, opts,is.num) {

  if(opts$outcome.class=="factor"){  #### THIS NEEDS TO BE UPDATED
    get.split <- bincartsplit(psi=psi, y=y, wt=wt, x.split=x.split, minsplit=minsplit, minbuck=minbuck, real.num=real.num, opts=opts,is.num=is.num)
  } else if (opts$outcome.class=="survival"){ #### THIS NEEDS TO BE UPDATED
    get.split <- survival.split(psi=psi, y=y, wt=wt, x.split=x.split,minsplit=minsplit,  minbuck=minbuck, real.num=real.num, opts=opts, is.num=is.num)
  } else {
    get.split <- cartsplit(psi=psi, y=y, wt=wt, x.split=x.split, minsplit=minsplit, minbuck=minbuck, real.num=real.num, opts=opts,is.num=is.num)
  }
  return(get.split)
}

cartsplit <- function(psi, y, wt, x.split, minsplit, minbuck, real.num, opts,is.num) {

  if(sum(is.na(x.split))>0){missing="yes"}else{missing="no"}
  
  ## psi tells which observations are in node that we are trying to split
  STOP <- 0
  x.s <- x.split[order(x.split)]
  wt.s <- wt[order(x.split)]
  y.s <- y[order(x.split)]
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
    y.s.k <- y.s.k[which(!is.na(x.s.k))]
    wt.s.k <- wt.s.k[which(!is.na(x.s.k))]
    x.s.k <- x.s.k[which(!is.na(x.s.k))] 
  }
 
  y.s.k <- y.s.k- sum(y.s.k * wt.s.k) / sum(wt.s.k)
  cant.split <- FALSE
  
  wtsum <- tapply(wt.s.k, x.s.k, sum)

  ## Added in the condition length(wt.s.k)!=0 since don't
  ## want to split a node with only one non-zero weight
  if(length(which(wtsum>0)) <= 1 | length(unique(y.s.k))==1 | length(y.s.k)<minsplit | length(which(wt.s.k!=0))<=1) {
    STOP <- 1
  } else if(abs(minsplit - max(table(x.s.k))) < minbuck &&
            ((length(x.s.k) - max(table(x.s.k))) < minbuck)) {
    			  STOP <- 1
  } else {  
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
    
    ## the denominator of goodness can sometimes be 0 which causes all
    ## goodness values to be NA. do we really need the denominator
    goodness <- (left.wt * lmean^2 + right.wt * rmean^2) / sum(wt.s.k * y.s.k^2)
	  if(sum(is.na(goodness))==length(goodness)){
         STOP <- 1
	  }else{
		  best.ind <- central.split(which.max(goodness),n)
      		max.g <- max(goodness, na.rm=TRUE)
      		in.set <- as.numeric(names(goodness)[best.ind])
  		if(missing=="no") {in.x <- (x.split <= (in.set + 1e-10))}
      		if(missing=="yes") {in.x <- (x.impute <= (in.set + 1e-10))}
      		## Need to consider number in the node where wt doesn't equal 0
      		num.in.set <- sum(in.x & psi==1 & wt >0)
      		num.not.in.set <- length(which(wt.s.k!=0)) - num.in.set
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
    
          ## Some of the temp values will be NA meaning some of the
          ## goodness values will be NA meaning that when we go through 
          ## the while loop, there is a chance that max.g could be NA as
          ## well and therefore I made this a condition in the while loop so
          ## the loop will keep repeating until we have a non-NA max.g value.
          while(is.na(max.g) | !(num.in.set >=  minbuck &&
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
            num.in.set <- sum(in.x & psi==1 & wt>0)
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
  if(STOP)  {
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
}


#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
bincartsplit <- function(psi, y, wt,  x.split, minsplit, minbuck, real.num, opts,is.num) {
  
  ## put everything in order by x so that we can find a split,
  ## ".s" stands for sorted
  STOP <- 0
  x.s <- x.split[order(x.split)]
  wt.s <- wt[order(x.split)]
  y.s <- y[order(x.split)]
  psi.s <- psi[order(x.split)]
  
  ## only look at the ones in this "OR bucket"
  x.s.k <- x.s[psi.s == 1] #.k indicates the OR bucket
  wt.s.k <- wt.s[psi.s == 1]
  y.s.k <- y.s[psi.s == 1]
  n.s.k <- length(y.s.k)
                                         # if there are missing values, impute the x to be x.impute using the mean or
                                        #  mode. then update x.s.k, y.s.k, wt.s.k to be the non-missing values since
                                        #  that is what will be used to create the splits. Because this is for a
                                        # categorical outcome, we average within each node for each outcome class.
  
  if(sum(is.na(x.split))>0){missing="yes"}else{missing="no"}
  if(missing=="yes"){
    
    x.original <- x.s.k
    x.impute<-x.split
    
    for (i in levels(y.s.k)){
      
                                        #warning scenario, but don't need warning message because it'll
                                        # already be outputted from update.missing
      if (length(which(!is.na(x.s.k[y.s.k==i])))==0){
        
        x.impute[which(is.na(x.split) & y==i)] <- ifelse(is.num==1,
                                          mean(x.s.k,na.rm=TRUE),as.numeric(names(which.max(table(x.s.k)))))
        
      }else{
        
        
        x.impute[which(is.na(x.split) & y==i)] <- ifelse(is.num==1,
                                          mean(x.s.k[y.s.k==i], na.rm=TRUE),
                                          as.numeric(names(which.max(table(x.s.k[y.s.k==i])))))
      }
    }
    y.s.k <- y.s.k[which(!is.na(x.s.k))]
    wt.s.k <- wt.s.k[which(!is.na(x.s.k))]
    x.s.k <- x.s.k[which(!is.na(x.s.k))]
    n.s.k <- length(y.s.k)
  }
  
  
  
  if(length(unique(y.s.k)) == 1 | length(unique(x.s.k))==1 | length(y.s.k)<minsplit)  {
    STOP <- 1
    
  } else{
    
    
    
    
    ux <- sort(unique(x.s.k))
    n=length(ux)
    
    parent.impurity <- gini.impurity(y.s.k, opts=opts)
    goodness <- rep(c(0), length(ux)-1)
    to.left<-0
    
    for(i in 1:(length(ux)-1)) {
      to.left<-to.left+sum(x.s.k==ux[i])
      goodness[i] <- parent.impurity -
        ((to.left/ n.s.k) * gini.impurity(y.s.k[1:to.left], opts=opts) +
         ((n.s.k - to.left) / n.s.k) * gini.impurity(y.s.k[(to.left + 1):n.s.k], opts=opts))
      
    }
    
    best.ind <- central.split(which.max(goodness),n)
    
    
    max.g <- max(goodness)
    
    if(max.g == 0) {
      STOP<-1
    
    }
    
    in.set <- ux[best.ind]
    if(missing=="no") {in.x <- (x.split <= (in.set + 1e-10))}
    if(missing=="yes") {in.x <- (x.impute <= (in.set + 1e-10))}
    num.in.set <- sum(in.x & psi==1)
    num.not.in.set <- length(x.s.k) - num.in.set
    
    
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
      
    } else{
      
      n.chg <- 2
      
      ## Given that best.ind are there enough on
      ## either side to split this region?
      while(!(num.in.set >= minbuck &&
              num.not.in.set >= minbuck)) {
        
        if(n.chg >= n) {
          STOP <- 1
          break
        }
        
        
        max.g <- goodness[order(goodness)][n - n.chg]
        n.chg <- n.chg + 1
        
        best.ind <- central.split(which(goodness == max.g),n)
        in.set <- ux[best.ind]
        
        if(missing=="no") {in.x <- (x.split <= (in.set + 1e-10))}
        if(missing=="yes") {in.x <- (x.impute <= (in.set + 1e-10))}
        
        num.in.set <- sum(in.x & psi==1)
        num.not.in.set <- length(x.s.k) - num.in.set
        
        
      } # end while statement
      
      if(missing=="no") {in.x <- (x.split <= (in.set + 1e-10))}
      if(missing=="yes") {in.x <- (x.impute <= (in.set + 1e-10))}
      
      num.in.set <- sum(in.x & psi==1)
      num.not.in.set <- length(x.s.k) - num.in.set
      
      new.lt <- as.numeric(in.x)
      new.gt <- 1 - new.lt
      new.lt <- ifelse(psi == 0, 0, new.lt)
      new.gt <- ifelse(psi == 0, 0, new.gt)
      val <- in.set + 1e-10
      cant.split <- FALSE
      
      
    } # end else
  } # end categorical x
  
  if(STOP) {   
    return(list(new.lt=psi, new.gt=rep(0, length(psi)),
                val=NA, cant.split=TRUE,
                max.g=(-1 * runif(1, 1000, 1010))))
  }
  
  
  return(list(new.lt=new.lt, new.gt=new.gt,
              val=val, cant.split=FALSE,
              max.g=max.g))
}
