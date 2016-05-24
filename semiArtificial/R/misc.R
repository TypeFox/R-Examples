maxValue <- .Machine$double.xmax # constant

# shift vector circularly for d elements, 
# shift left with negative d, shift right with positive d
shift <- function(x, d=1) {
  right <- TRUE # positive d means shift right
  if (d < 0) {
    right <- FALSE
    d <- -d
  }
  d <- d %% length(x)
  if (d == 0)
    return(x)
  if (right) {  
    start <-x[1:(length(x)-d)]
    end <- x[(length(x)-d+1):length(x)]
  }
  else {
    start <-x[1:d]
    end <- x[(d+1):length(x)]
  }
  c(end,start)    
}

resample <- function(x, ...) x[sample.int(length(x), ...)]

sqr <- function(x) {
  x*x
}

#generate k-fold cross validation of n instances
cvGen<-function(n,k) {
  v <- 1:k
  vec <- array(1:k,dim=n)
  sample(vec, size=n)
}

# generate stratified k-fold cross validation partition based on classes in classVal
cvGenStratified<-function(classVal,k) {
  classVal<-factor(classVal)
  levs = factor(levels(classVal), levels=levels(classVal))
  classFreq <- table(classVal)
  noClasses <- length(levs)
  n <- length(classVal)
  srt <- order(classVal)
  vec <- array(1:k,dim=n)
  cv <- array(0,dim=n)
  cv[srt]<-vec
  for (i in 1:noClasses) 
    cv[classVal==levs[i]] <- sample(cv[classVal==levs[i]], size=classFreq[i], replace=F)
  cv
}

# collect instances with the same position in different sublists of lst
gatherFromList<-function(lst){
  m <-list()
  
  for (j in 1:length(lst[[1]])) {
    m[[j]]<-vector(mode="numeric",length=length(lst))
    names(m[[j]]) <- names(lst)
  }
  names(m)<-names(lst[[1]])
  for (i in 1:length(lst)){
    for (j in 1:length(lst[[i]])){
      if (is.null(dim(lst[[i]][[j]])))
        m[[j]][i] <- lst[[i]][[j]]
    }
  }  
  m
}

# collect instances with the same position in different components of lst
dataframeFromList<-function(lst){
  df <- data.frame(x = rep(NA, length(lst)))
    
  for (j in 1:length(names(lst[[1]]))) {
    col <- c() 
    for (i in 1:length(lst)){
      col[i] <-lst[[i]][[j]]  
    }
    df[[j]] <- col
  }  
  #names(df) <- names(lst[[1]])
  df
}



# generator for monmk 1,2,3 data sets
monkGen<-function(noInst, problem=1, pYes=0.5, classNoise=0) {
  missingYes<-as.integer(noInst*pYes)
  missingNo<-noInst-missingYes
  m<-data.frame()
  while (missingYes>0 || missingNo > 0) {
    D1 <- as.integer(runif(noInst,0,3))
    D2 <- as.integer(runif(noInst,0,3))
    D3 <- as.integer(runif(noInst,0,2))
    D4 <- as.integer(runif(noInst,0,3))
    D5 <- as.integer(runif(noInst,0,4))
    D6 <- as.integer(runif(noInst,0,2))
    f<-numeric(noInst)
    if (problem==1) {
      f[D1==D2 | D5==0] <- 1
    }
    else if (problem==2) {
      f[D1==0] <- f[D1==0]+1
      f[D2==0] <- f[D2==0]+1
      f[D3==0] <- f[D3==0]+1
      f[D4==0] <- f[D4==0]+1
      f[D5==0] <- f[D5==0]+1
      f[D6==0] <- f[D6==0]+1
      f[f!=2] <- 0
      f[f==2] <- 1    
    }
    else if (problem==3) {
      f[(D5==2 & D4==0) | (D5!=3 & D2!=2)] <- 1
    }
    else stop("Invalid number of monk problem, should be 1, 2, or 3")
    #class noise
    if (classNoise > 0) {
      ns<-runif(noInst,0,1)
      f[ns<classNoise] <- abs(f[ns<classNoise]-1) # revert
    }
    t <- data.frame(
      head_shape=factor(D1),
      body_shape=factor(D2),
      is_smiling=factor(D3),
      holding=factor(D4),
      jacket_color=factor(D5),
      has_tie=factor(D6),
      class=factor(f)
    )
    cl1=t[t$class==1,]
    cl0=t[t$class==0,]
    if (nrow(cl1)>=missingYes) {
      m<-rbind(m,cl1[seq(1,length.out=missingYes),])
      missingYes<-0
    }
    else {
      m<-rbind(m,cl1)
      missingYes<-missingYes-nrow(cl1)
    }
    if (nrow(cl0)>=missingNo) {
      m<-rbind(m,cl0[seq(1,length.out=missingNo),])
      missingNo<-0
    }
    else {
      m<-rbind(m,cl0)
      missingNo<-missingNo-nrow(cl0)
    }
  }
  
  levels(m$head_shape)<-c("round","square","octagon")
  levels(m$body_shape)<-c("round","square","octagon")
  levels(m$is_smiling)<-c("yes","no")
  levels(m$holding)<-c("sword","baloon","flag")
  levels(m$jacket_color)<-c("red","yellow","green","blue")
  levels(m$has_tie)<-c("yes","no")  
  levels(m$class)<-c("no","yes")  
  
  monk<-m[sample(noInst),]
  monk
}

# the data.frame set is converted to a form such that all the attributes have values between 0 and 1
# this is useful in visualization
varNormalization<-function(md, set){
  #d-discrete, a-attribure, n-names
  column<-length(set[1,]);
  n<-length(set[,1]);
  colPos<-matrix(FALSE, column, column);
  dan<-md$discAttrNames;
  nd<-length(dan);
  ian<-0;
  if(nd>0){
    int<-vector("numeric",nd);
    for(ian in 1:nd)
    {
      search<-dan[ian];
      colPos[ian,]<-unlist(lapply(attr(set, 'names'),function(x){x==search}));
      int[ian]<-1/(length(levels(set[, colPos[ian, ]]))-1);
    }
  }
  nan<-md$numAttrNames;
  nn<-length(nan);
  if(nn > 0){
    offset<-ian;
    mi<-vector("numeric", nn);
    sigma<-vector("numeric", nn);
    maxnorm<-vector("numeric", nn);
    moveup<-vector("numeric", nn);
    for(ian in 1:nn){
      search<-nan[ian];
      tmp<-unlist(lapply(attr(set, 'names'),function(x){x==search}));
      colPos[ian+offset,]<-tmp
      mi[ian]<-as.numeric(mean(set[,tmp]));
      sigma[ian]<-as.numeric(sd(set[,tmp]));
      allcurval<-(set[, tmp]-mi[ian])/sigma[ian];
      moveup[ian]<-min(c(0,allcurval));
      maxnorm[ian]<-max(allcurval-moveup[ian])
    }
  }
  out<-NULL;
  classV<- attr(md$terms, "variables")[[2]];
  classV<-unlist(lapply(attr(set, 'names'),function(x){x!=classV}));
  for(ex in 1:n){
    pos<-vector("numeric", column);
    if(nd>0){
      for(da in 1:nd) {
        lev<-levels(set[, colPos[da, ]]);
        val<-set[ex, colPos[da, ]];
        if(is.na(val)){
          index<-1;
        }
        else{
          index<-c(1:length(lev))[lev==val];
        }
        pos[colPos[da, ]]<-(index-1)*int[da];
      }
    }
    if(nn>0){
      for(na in 1:nn){
        normal<-((set[ex, colPos[na+offset, ]]-mi[na])/sigma[na]);
        val<-(normal-moveup[na])/maxnorm[na];
        pos[colPos[na+offset, ]]<-val;
      }
    }
    out<-c(out,pos[classV]);
  }
  out<-matrix(out, n, column-1, TRUE);
}

plotRFNorm<-function(point, cluster, somnames, lOffset,myHoriz=FALSE, myAxes=FALSE)
{
  noVar<-dim(point)[2];
  ylim<-c(min(point), max(point)+lOffset)
  xlim<-c(1, noVar)
  if(is.logical(myAxes)){
    axesShow<-TRUE;
  }
  else{
    axesShow<-FALSE;
  }
  plot(1, 1, xlim=xlim, ylim=ylim, type="n", ann=FALSE, frame=TRUE, axes=axesShow);
  if(!is.logical(myAxes) && length(myAxes) > 0){
    axis(2);
    axis(1, at=1:noVar, labels=myAxes);
  }
  tmpCluster <- cluster;
  clusterLevelNames <- list();
  clusterLevels <- 0;
  i <- 1;
  while(length(tmpCluster) > 0 && i < 13){
    clusterLevels[i] <- i;
    clusterLevelNames[i]<-tmpCluster[1]
    cluster[cluster==tmpCluster[1]]<-i;
    tmpCluster <- tmpCluster[tmpCluster!=tmpCluster[1]];
    i <- i+1;
  }
  clusterLevels <- clusterLevels[clusterLevels>0];
  myPch <-0;
  myColor <-0;
  myCount<-7
  nexamples<-length(point[,1]);
  for(value in 1:nexamples){
    #mod
    myPch[cluster[value]]<-floor(cluster[value]/myCount);
    myColor[cluster[value]]<- 1+cluster[value] - myCount*floor(cluster[value]/myCount);
    points(point[value,], col=myColor[cluster[value]], pch=myPch[cluster[value]]);
    lines(point[value,], col=myColor[cluster[value]], pch=myPch[cluster[value]]);}
  legend(xlim[1], ylim[2], somnames, cex=0.8, col=myColor, pch=myPch, horiz=myHoriz);
}

eps <- 1e-6 

#computes cdf from collection of probabilities summing to 1
probs2cdf<-function(probs) {
  cdf <- probs
  i <- 2
  while (i <= length(probs)){
    cdf[i] <- cdf[i-1] + cdf[i]
    i <- i+1
  }
  if (abs(cdf[length(cdf)] - 1.0)  > eps) 
    warning("Sum of probability distribution is not equal to 1")
  cdf[length(cdf)] <- 1.0
  i <- length(cdf) - 1
  while (cdf[i] > 1.0) {
    cdf[i] <- 1.0
    i <- i -1
  }
  cdf
}

# finds quantile of a factor given by probability distribution
quant<-function(p, probs) {
  cdf <- probs2cdf(probs)
  val <- findInterval(p, cdf, rightmost.closed=T)
  return(val+1)
}

probs2str<-function(x, digits=2) {
  s <- "["
  for (i in seq(along.with=x)){
    s <- paste(s,format(x[i],digits=digits),sep="")
    if (i < length(x))
      s<-paste(s,", ",sep="")
  }
  paste(s,"]",sep="")
}

factors2str<-function(x) {
  s <- "["
  for (i in seq(along.with=x)){
    s <- paste(s,format(x[i]),sep="")
    if (i < length(x))
      s<-paste(s,", ",sep="")
  }
  paste(s,"]",sep="")
}


intFromProb <- function(probs, n) {
  fairN <- probs * n
  roundN <- round(fairN)
  sumN <- sum(roundN,na.rm=T)
  while (sumN > n) {
    mostUnfair <- which.is.max(roundN-fairN)
    roundN[mostUnfair] <- roundN[mostUnfair] -1
    sumN <- sumN - 1
  }
  while (sumN < n) {
    mostUnfair <- which.is.max(fairN-roundN)
    roundN[mostUnfair] <- roundN[mostUnfair] + 1
    sumN <- sumN + 1
  }
  roundN   
}

#imputation with median value
medianImpute <- function(data) {
  for (i in 1:ncol(data)) {
    nas <- is.na(data[[i]])
    if (any(nas)) {
      if (is(data[[i]],"factor"))
        data[which(nas), i] <- factorMode(data[[i]])
      else
       data[which(nas), i] <- median(data[[i]],na.rm=T)
    }
  }
  data
}

# modified median - always returns actual value
mmedian<- function(x) {
  sx <- sort(x)
  sx[length(sx)/2]
}

# mode - the most frequent value of a factor
factorMode <- function(x) {
  freq <- table(x)
  f<-names(freq)[which.is.max(freq)]
  factor(f, levels=levels(x))
}

normalize01 <- function(data) {
  for (i in 1:ncol(data)) {
    cmin <- min(data[,i], na.rm=TRUE)
    cmax <- max(data[,i], na.rm=TRUE)
    if (cmin < cmax)
       data[,i] <- (data[,i]-cmin)/(cmax-cmin)
    else data[,i] <- 0
  }
  data
}
  
# computes Hellinger distance between two discrete distributions given in the form of probability tables
hellinger <- function(dist1, dist2) {
  dist1 <- sqrt(dist1)
  dist2 <- sqrt(dist2)
  dd <- (dist1-dist2)^2
  h <- sqrt(sum(dd))/sqrt(2)
  return(h)
}

# a simple toy data generator, producing separated Gaussians, with two or three class values and one to three attributes
toyGauss <- function(noInst=100, noAttr=2, noClassValues=3){
	noData <- noInst
    if (! noClassValues %in% c(2,3))
      stop("Argument noClassValues shall be 2 or 3.")
	if (! noAttr %in% 1:3)
		stop("Argument noAttr shall be 1, 2, or 3.")
	
	mu1 <- rep(-5,noAttr)
	sigma1 <- diag(nrow=noAttr)
	attrs1 <- mvrnorm(n=noData,mu=mu1,Sigma=sigma1)
	data1<-cbind(attrs1,rep(1,noData)) # attach class value 1
	
	mu2 <- rep(5,noAttr)
	sigma2 <- diag(nrow=noAttr)
	attrs2 <- mvrnorm(n=noData,mu=mu2,Sigma=sigma2)
	data2<-cbind(attrs2,rep(2,noData)) # attach class value 2
	
	mu3 <- rep(0,noAttr)
	sigma3 <- diag(nrow=noAttr)
	attrs3 <- mvrnorm(n=noData,mu=mu3,Sigma=sigma3)
	data3<-cbind(attrs3,rep(3,noData)) # attach class value 3
	
	if (noClassValues==3)
		data <- rbind(data1,data2,data3) # merge 3
	else if (noClassValues==2)
		data <- rbind(data1,data2) # merge 2
	
	data <-as.data.frame(data)
	#form attribute names
	attrNames <- c("class")
	attrNames <- c(paste("A",1:noAttr,sep=""),attrNames)
	names(data) <- attrNames
	data$"class" <- as.factor(data$"class")
	data
}
