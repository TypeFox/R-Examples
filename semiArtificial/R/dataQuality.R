# functions comparing the data sets and measuring quality of generated data


# returns different statistics comparing the similarity of data1 and data2
# some statistics are applied to discrete and other too numeric columns only
# dropDiscrete contains a list of discrete features which shall be left out of comparison i.e., response variable
dataSimilarity <- function(data1, data2, dropDiscrete=NA) {
  
  if(ncol(data1)!=ncol(data2))
    stop("Only data with equal number of columns can be compared.") 
  
  # identify discrete columns
  discrete <- c()
  noVal <- c()
  for (i in 1:ncol(data1)) {
    if (is.factor(data1[[i]]) != is.factor(data2[[i]]))
      stop("Only data with equal types of columns can be compared.", names(data1)[i], " ", names(data2)[i])
    discrete[i] <- is.factor(data1[[i]])
    if (discrete[i]) {
      if (! all(levels(data1[[i]]) == levels(data2[[i]])))
        stop("Only data with equal types of columns can be compared.", levels(data1[[i]]), " ", levels(data2[[i]]))
      noVal[i] <- length(levels(data1[[i]]))
    }
    else
      noVal[i] <- NA
  }
  
  # extract numeric columns
  data1Num <- data1[,!discrete,drop=F]
  data2Num <- data2[,!discrete,drop=F]
  
  #prepare result matrices
  s1<- matrix(0,nrow=7,ncol=ncol(data1Num))
  s2<- matrix(0,nrow=7,ncol=ncol(data1Num))
  s1Norm<- matrix(0,nrow=7,ncol=ncol(data1Num))
  s2Norm<- matrix(0,nrow=7,ncol=ncol(data1Num))
  ks <- c() # store p-values of Kolmogorov-Smirnov test on drawing from the same distribution, only for numeric attributes
  
  # numeric columns
  if (ncol(data1Num) > 0) {
    
    s1[1,] <- apply(data1Num, 2, mean, na.rm=TRUE)
    s1[2,] <- apply(data1Num, 2, sd, na.rm=TRUE)
    s1[3,] <- apply(data1Num, 2, skewness, na.rm=TRUE)
    s1[4,] <- apply(data1Num, 2, kurtosis, na.rm=TRUE)
    s1[5,] <- apply(data1Num, 2, mc, na.rm=TRUE) # medcouple MC
    med1 <- apply(data1Num, 2, median, na.rm=TRUE)
    for (i in 1:ncol(data1Num)) {
      s1[6,i] <-  mc(data1Num[data1Num[,i] < med1[i], i], na.rm=TRUE) # left MC, LMC
      s1[7,i] <-  mc(data1Num[data1Num[,i] > med1[i], i], na.rm=TRUE) # right MC, RMC
    }
    
    s2[1,] <- apply(data2Num, 2, mean, na.rm=TRUE)
    s2[2,] <- apply(data2Num, 2, sd, na.rm=TRUE)
    s2[3,] <- apply(data2Num, 2, skewness, na.rm=TRUE)
    s2[4,] <- apply(data2Num, 2, kurtosis, na.rm=TRUE)
    s2[5,] <- apply(data2Num, 2, mc, na.rm=TRUE) # medcouple MC
    med2 <- apply(data2Num, 2, median, na.rm=TRUE)
    for (i in 1:ncol(data2Num)) {
      s2[6,i] <-  mc(data2Num[data2Num[,i] < med2[i], i], na.rm=TRUE) # left MC, LMC
      s2[7,i] <-  mc(data2Num[data2Num[,i] > med2[i], i], na.rm=TRUE) # right MC, RMC
    }
    # normalize data to [0,1] 
    data1NumNorm <- normalize01(data1Num)
    data2NumNorm <- normalize01(data2Num)
    
    s1Norm[1,] <- apply(data1NumNorm, 2, mean, na.rm=TRUE)
    s1Norm[2,] <- apply(data1NumNorm, 2, sd, na.rm=TRUE)
    s1Norm[3,] <- apply(data1NumNorm, 2, skewness, na.rm=TRUE)
    s1Norm[4,] <- apply(data1NumNorm, 2, kurtosis, na.rm=TRUE)
    s1Norm[5,] <- apply(data1NumNorm, 2, mc, na.rm=TRUE) # medcouple MC
    med1n <- apply(data1NumNorm, 2, median, na.rm=TRUE)
    for (i in 1:ncol(data1NumNorm)) {
      s1Norm[6,i] <-  mc(data1NumNorm[data1NumNorm[,i] < med1n[i], i], na.rm=TRUE) # left MC, LMC
      s1Norm[7,i] <-  mc(data1NumNorm[data1NumNorm[,i] > med1n[i], i], na.rm=TRUE) # right MC, RMC
    }
    
    s2Norm[1,] <- apply(data2NumNorm, 2, mean, na.rm=TRUE)
    s2Norm[2,] <- apply(data2NumNorm, 2, sd, na.rm=TRUE)
    s2Norm[3,] <- apply(data2NumNorm, 2, skewness, na.rm=TRUE)
    s2Norm[4,] <- apply(data2NumNorm, 2, kurtosis, na.rm=TRUE)
    s2Norm[5,] <- apply(data2NumNorm, 2, mc, na.rm=TRUE) # medcuple MC
    med2n <- apply(data2NumNorm, 2, median, na.rm=TRUE)
    for (i in 1:ncol(data2NumNorm)) {
      s2Norm[6,i] <-  mc(data2NumNorm[data2NumNorm[,i] < med2n[i], i], na.rm=TRUE) # left MC, LMC
      s2Norm[7,i] <-  mc(data2NumNorm[data2NumNorm[,i] > med2n[i], i], na.rm=TRUE) # right MC, RMC
    }
    
    warn.save <- getOption("warn")
    options(warn=-1) # switch off warnings about ties  
    for (i in 1:ncol(data1Num)) { # KS test 
      tst <- ks.test(data1Num[,i], data2Num[,i])
      ks[i] <- tst$p.value
    }
    names(ks) <- names(data1Num)
    options(warn=warn.save) # restore warning level
    
  }
  # difference between normalized statistics
  ds <- s1Norm - s2Norm
  
  colnames(s1) <- colnames(s2) <- colnames(s1Norm) <- colnames(s2Norm) <- colnames(ds) <- names(data1Num)
  rownames(s1) <- rownames(s2) <- rownames(s1Norm)<- rownames(s2Norm) <- rownames(ds) <- c("mean","stdev","skewness","kurtosis","medcouple","LMC","RMC")
  
  # nominal attributes: compare frequencies and Hellinger distance
  f1 <- list()
  f2 <- list()
  df <- list()
  hel <- c() # store Hellinger distance between two discrete distributions
  discreteIdx<-which(discrete)
  if (! is.na(dropDiscrete)) {
    dropIdx <- match(dropDiscrete, names(data1))
    discreteIdx <- discreteIdx[-match(dropIdx, discreteIdx, nomatch=length(discreteIdx)+1)]
  }
  for (i in seq(along=discreteIdx)) {
    f1[[i]] <- table(data1[,discreteIdx[i]])/nrow(data1)
    names(f1)[i] <-names(data1)[discreteIdx[i]]
    f2[[i]] <- table(data2[,discreteIdx[i]])/nrow(data2)
    names(f2)[i] <-names(data2)[discreteIdx[i]]
    df[[i]] <- f1[[i]] - f2[[i]]
    names(df)[i] <-names(data1)[discreteIdx[i]]
    hel[i] <- hellinger(f1[[i]], f2[[i]])
  }
  if (length(hel) > 0) 
    names(hel)<- names(data1)[discreteIdx]
  
  # return values
  list(equalInstances=noEqualRows(data1,data2), stats1num=s1, stats2num=s2, stats1numNorm=s1Norm, stats2numNorm=s2Norm,
       KSpvalue=ks, freq1=f1, freq2=f2, dfreq=df, dstatsNorm=ds, hellingerDist=hel)
}


# find a number of the same instances in data1 and data2 (up to the given tolerance)
# use either a double loop or construct a distance matrix depending on data dimensions
equalInst<-function(data1, data2) {
  if (nrow(data1) > ncol(data1))
    return(equalInstRow(data1, data2))
  else
    return(equalInstDaisy(data1, data2))
}


equalInstRow<-function(data1, data2) {
  eq <- 0
  for (i in 1:nrow(data1)) 
    for (j in 1:nrow(data2))    
      eq <- eq + as.integer(almostEqual(data1[i,], data2[j,],tolerance=1e-5))
  eq
}

# find a number of the same instances in data1 and data2 (up to the given tolerance)
# use the Gower's metric and daisy distance matrix
equalInstDaisy <- function(data1, data2, tolerance=1e-5) {
  dta <- rbind(data1, data2)
  sim <- as.matrix(daisy(dta, metric="gower"))
  eq <- 0
  for (i in (nrow(data1)+1):nrow(sim)) {
    if (length(which(sim[i,]< tolerance)) >1)
      eq <- eq+1
  }
  eq
}

# for each instance in data2 find a closest instance in data1 and 
# compute a mean of these minimal distances
meanMinDistance<-function(data1, data2) {
  dist <- vector(mode="numeric",length=nrow(data2))
  for (i in 1:nrow(data2)) {
    row2 <- as.numeric(data2[i,])
    dist[i] <- min(apply(data1, 1, diff, row2))
  }
  
  list(mean=mean(dist),min=min(dist),dist=dist)
}

# sum of squared distance between two instances
diff<-function(inst1, inst2){
  sqrt(sum((inst1-inst2)^2,na.rm=T))
}

# are two instances almost equal 
almostEqual <- function(inst1, inst2, tolerance=1e-5) {
  all(abs(as.numeric(inst1)-as.numeric(inst2)) < tolerance, na.rm=TRUE)
}

# compare adjusted Rand index and/or Fowlkes-Mallows index of two cluserings obtained by:
# 1. using data1 to get clusters, assign data2 to these clusters
# 2. using data2 to get clusters, assign data1 to these clusters
# The method uses pam clustering (partitioning around medoids) with 
# optimal k computed for each data1 and data2 separately with average silhouette width method (using pamk)
# 
# ARI is computed on merged data1 and data2 instances 
ariCompare <- function(data1, data2, ARI=TRUE, FM=TRUE) {
  
  # check if data is of compatibile types
  if(ncol(data1)!=ncol(data2))
    stop("Only data with equal number of columns can be compared.")
  discrete <- c()
  noVal <- c()
  for (i in 1:ncol(data1)) {
    if (is.factor(data1[[i]]) != is.factor(data2[[i]]))
      stop("Only data with equal types of columns can be compared.", names(data1)[i], " ", names(data2)[i])
    discrete[i] <- is.factor(data1[[i]])
    if (discrete[i]) {
      if (! all(levels(data1[[i]]) == levels(data2[[i]])))
        stop("Only data with equal types of columns can be compared.", levels(data1[[i]]), " ", levels(data2[[i]]))
      noVal[i] <- length(levels(data1[[i]]))
    }
    else
      noVal[i] <- NA
  }
  
  #distance matrix
  dist1 <- daisy(data1, metric="gower")
  dist2 <- daisy(data2, metric="gower")
  
  pamkBest1 <- pamk(dist1) # number of clusters estimated by optimum average silhouette width
  pamkBest2 <- pamk(dist2) # number of clusters estimated by optimum average silhouette width
  
  clu1 <- pam(dist1, k=max(2,pamkBest1$nc)) #get the clustering
  clu2 <- pam(dist2, k=max(2,pamkBest2$nc)) #get the clustering
  
  medoids1 <- data1[clu1$medoids,] 
  medoids2 <- data2[clu2$medoids,] 
  
  # compute distances of data2 to medoids1 and data1 to medoids2
  dta1 <- rbind(medoids1,data2)
  dta2 <- rbind(medoids2,data1)
  dst1 <- as.matrix(daisy(dta1,metric="gower"))
  dst2 <- as.matrix(daisy(dta2,metric="gower"))
  # find nearest medoid for each data2 instance, this is its cluster index in 1st clustering
  closestMed1 <-apply(dst1[-c(1:nrow(medoids1)),1:nrow(medoids1)],1,which.min)
  # find nearest medoid for each data1 instance, this is its cluster index in 2nd clustering
  closestMed2 <-apply(dst2[-c(1:nrow(medoids2)),1:nrow(medoids2)],1,which.min)
  
  # put the cluster assignments together
  clustering1 <- c(clu1$clustering, closestMed1)
  clustering2 <- c(closestMed2, clu2$clustering)
  
  if (ARI)
    ari <- adjustedRandIndex(clustering1, clustering2)
  else ari <- NULL
  if (FM)
    fm <- FM_index(clustering1, clustering2)
  else fm <- NULL
  
  list(ARI=ari, FM=fm)
  
}

# compare classification performance on two data sets of the same problem
# by splitting the each data set into two parts and do a 2-fold stratified cross-validation
performanceCompare <- function(data1, data2, formula, model="rf", stat=NULL, ...) {
  
  informula <- as.formula(formula)
  response <- all.vars(informula)[1]
  
  # check if data is of compatibile types
  if(ncol(data1)!=ncol(data2))
    stop("Only data with equal number of columns can be compared.")
  discrete <- c()
  noVal <- c()
  for (i in 1:ncol(data1)) {
    if (is.factor(data1[[i]]) != is.factor(data2[[i]]))
      stop("Only data with equal types of columns can be compared.", names(data1)[i], " ", names(data2)[i])
    discrete[i] <- is.factor(data1[[i]])
    if (discrete[i]) {
      if (! all(levels(data1[[i]]) == levels(data2[[i]])))
        stop("Only data with equal types of columns can be compared.", levels(data1[[i]]), " ", levels(data2[[i]]))
      noVal[i] <- length(levels(data1[[i]]))
    }
    else
      noVal[i] <- NA
  }
  if (is.factor(data1[[response]]))
    problemType <- "classification"
  else problemType <- "regression"
  if (is.null(stat))
    if (problemType == "classification")
        stat <- "accuracy"
    else stat <- "RMSE"
  
  # shuffle the data to allow several different runs
  d1 <- data1[sample(nrow(data1),nrow(data1)),]
  d2 <- data2[sample(nrow(data2),nrow(data2)),]
  
  if (problemType == "classification") {
    # split into halves based on startified cross-validation
    cv1 <- cvGenStratified(d1[[response]], k=2)
    d1a <- d1[cv1==1,]
    d1b <- d1[cv1==2,]
    
    cv2 <- cvGenStratified(d2[[response]], k=2)
    d2a <- d2[cv2==1,]
    d2b <- d2[cv2==2,]
    
    # build models
    m1a <- CoreModel(informula, d1a, model=model, ...)
    m1b <- CoreModel(informula, d1b, model=model, ...)
    
    m2a <- CoreModel(informula, d2a, model=model, ...)
    m2b <- CoreModel(informula, d2b, model=model, ...)
    
    # test performance
    # model m1a
    pred.m1ad1b <- predict(m1a, d1b)
    acc.m1ad1b <- modelEval(m1a, d1b[[response]], pred.m1ad1b$class, pred.m1ad1b$probabilities)  
    
    pred.m1ad2 <- predict(m1a, d2)
    acc.m1ad2 <- modelEval(m1a, d2[[response]], pred.m1ad2$class, pred.m1ad2$probabilities)  
    
    #model m1b
    pred.m1bd1a <- predict(m1b, d1a)
    acc.m1bd1a <- modelEval(m1b, d1a[[response]], pred.m1bd1a$class, pred.m1bd1a$probabilities)  
    
    pred.m1bd2 <- predict(m1b, d2)
    acc.m1bd2 <- modelEval(m1b, d2[[response]], pred.m1bd2$class, pred.m1bd2$probabilities)  
    
    # model m2a
    pred.m2ad2b <- predict(m2a, d2b)
    acc.m2ad2b <- modelEval(m2a, d2b[[response]], pred.m2ad2b$class, pred.m2ad2b$probabilities)  
    
    pred.m2ad1 <- predict(m2a, d1)
    acc.m2ad1 <- modelEval(m2a, d1[[response]], pred.m2ad1$class, pred.m2ad1$probabilities)  
    
    #model m2b
    pred.m2bd2a <- predict(m2b, d2a)
    acc.m2bd2a <- modelEval(m2b, d2a[[response]], pred.m2bd2a$class, pred.m2bd2a$probabilities)  
    
    pred.m2bd1 <- predict(m2b, d1)
    acc.m2bd1 <- modelEval(m2b, d1[[response]], pred.m2bd1$class, pred.m2bd1$probabilities)  
    
    # destroy models
    destroyModels(c(m1a,m1b,m2a,m2b))
  }
  else { # regression
    # split into halves based on cross-validation
    cv1 <- cvGen(n=nrow(d1), k=2)
    d1a <- d1[cv1==1,]
    d1b <- d1[cv1==2,]
    
    cv2 <- cvGen(n=nrow(d2), k=2)
    d2a <- d2[cv2==1,]
    d2b <- d2[cv2==2,]
    
    # build models
    if (model == "regTree") {
      m1a <- CoreModel(informula, d1a, model=model, ...)
      m1b <- CoreModel(informula, d1b, model=model, ...)
      
      m2a <- CoreModel(informula, d2a, model=model, ...)
      m2b <- CoreModel(informula, d2b, model=model, ...)
    }
    else if (model %in% c("rf","bagging")) {
      m1a <- do.call(model, list(informula, d1a, ...))
      m1b <- do.call(model, list(informula, d1b, ...))
      
      m2a <- do.call(model, list(informula, d2a, ...))
      m2b <- do.call(model, list(informula, d2b, ...))
    }
    # test performance
    # model m1a
    pred.m1ad1b <- predict(m1a, d1b)
    acc.m1ad1b <- modelEval(NULL, d1b[[response]], pred.m1ad1b, avgTrainPrediction=mean(d1a[[response]]))  
    
    pred.m1ad2 <- predict(m1a, d2)
    acc.m1ad2 <- modelEval(NULL, d2[[response]], pred.m1ad2, avgTrainPrediction=mean(d1a[[response]]))  
    
    #model m1b
    pred.m1bd1a <- predict(m1b, d1a)
    acc.m1bd1a <- modelEval(NULL, d1a[[response]], pred.m1bd1a, avgTrainPrediction=mean(d1b[[response]]))  
    
    pred.m1bd2 <- predict(m1b, d2)
    acc.m1bd2 <- modelEval(NULL, d2[[response]], pred.m1bd2, avgTrainPrediction=mean(d1b[[response]]))  
    
    # model m2a
    pred.m2ad2b <- predict(m2a, d2b)
    acc.m2ad2b <- modelEval(NULL, d2b[[response]], pred.m2ad2b, avgTrainPrediction=mean(d2a[[response]]))  
    
    pred.m2ad1 <- predict(m2a, d1)
    acc.m2ad1 <- modelEval(NULL, d1[[response]], pred.m2ad1, avgTrainPrediction=mean(d2a[[response]]))  
    
    #model m2b
    pred.m2bd2a <- predict(m2b, d2a)
    acc.m2bd2a <- modelEval(NULL, d2a[[response]], pred.m2bd2a, avgTrainPrediction=mean(d2b[[response]]))  
    
    pred.m2bd1 <- predict(m2b, d1)
    acc.m2bd1 <- modelEval(NULL, d1[[response]], pred.m2bd1, avgTrainPrediction=mean(d2b[[response]]))  
    
    # destroy CoreModel models
    if (model == "regTree")
		destroyModels(c(m1a,m1b,m2a,m2b))
	else{
		m1a <- m1b <- m2a <- m2b <- NULL
	}  
  }
  # averages for 2-fold cross-validation
  acc.m1d1 <- (acc.m1ad1b[[stat]] + acc.m1bd1a[[stat]])/2.0
  acc.m1d2 <- (acc.m1ad2[[stat]] + acc.m1bd2[[stat]])/2.0
  acc.m2d1 <- (acc.m2ad1[[stat]] + acc.m2bd1[[stat]])/2.0
  acc.m2d2 <- (acc.m2ad2b[[stat]] + acc.m2bd2a[[stat]])/2.0
  diff.m1 <- acc.m1d1 -acc.m1d2
  diff.m2 <- acc.m2d2 - acc.m2d1

  list(diff.m1=diff.m1, diff.m2=diff.m2, perf.m1d1=acc.m1d1, perf.m1d2=acc.m1d2, perf.m2d1=acc.m2d1, perf.m2d2=acc.m2d2 )  
}

# compare classification performance on two data sets of the same problem
performanceCompareSingle <- function(data1, data2, formula, model="rf", stat="accuracy", ...) {
  
  informula <- as.formula(formula)
  trms <- terms(as.formula(informula),data=data1)
  response <- all.names(attr(trms,"variables"))[1+attr(trms,"response")]
  
  # check if data is of compatibile types
  if(ncol(data1)!=ncol(data2))
    stop("Only data with equal number of columns can be compared.")
  discrete <- c()
  noVal <- c()
  for (i in 1:ncol(data1)) {
    if (is.factor(data1[[i]]) != is.factor(data2[[i]]))
      stop("Only data with equal types of columns can be compared.", names(data1)[i], " ", names(data2)[i])
    discrete[i] <- is.factor(data1[[i]])
    if (discrete[i]) {
      if (! all(levels(data1[[i]]) == levels(data2[[i]])))
        stop("Only data with equal types of columns can be compared.", levels(data1[[i]]), " ", levels(data2[[i]]))
      noVal[i] <- length(levels(data1[[i]]))
    }
    else
      noVal[i] <- NA
  }
  
  # shuffle the data
  d1 <- data1[sample(nrow(data1),nrow(data1)),]
  d2 <- data2[sample(nrow(data2),nrow(data2)),]
  
  
  # build models
  m1 <- CoreModel(informula, d1, model=model, ...)
  m2 <- CoreModel(informula, d2, model=model, ...)
  
  # test performance
  # model m1
  pred.m1d2 <- predict(m1, d2)
  acc.m1d2 <- modelEval(m1, d2[[response]], pred.m1d2$class, pred.m1d2$probabilities)  
  
  # model m2a
  pred.m2d1 <- predict(m2, d1)
  acc.m2d1 <- modelEval(m2, d1[[response]], pred.m2d1$class, pred.m2d1$probabilities)  
  
  # destroy models
  destroyModels(c(m1,m2))
  
  # difference in performance
  diff <- acc.m1d2[[stat]] - acc.m2d1[[stat]]
  
  list(m1=acc.m1d2[[stat]], m2=acc.m2d1[[stat]], diff=diff)
}


