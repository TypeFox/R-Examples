#####
## Collection of internal functions used
##
#####
.onAttach <- function(...)
   {
   packageStartupMessage(paste(" **********   PK Version",
packageDescription("PK")$Version), "********** \n")
  # cat(paste("For more on PK see the vignette.\n"))
   packageStartupMessage("Type PKNews() to see new features/changes/bug fixes.\n")
   if(as.numeric(R.Version()$major) <= 2 & as.numeric(R.Version()$minor) < 10) {packageStartupMessage(paste("The help functions for this package might not work properly. Please upgrade R to Version 2.10 or above to fix this problem.\n"))}
   }


# function to calculate weights for linear trapezoidal rule
.weight <- function(time){
  if(length(time) == 0){return(0)}
  time <- unique(sort(unlist(time)))
  J <- length(time)
  i <- c(2:(J-1))
  w <- rep(NA,J)
  w[1] <- (time[2] - time[1])/2
  w[i] <- (time[i+1] - time[i-1])/2
  w[J] <- (time[J] - time[J-1])/2
  return(w)
}

# Fieller confidence interval for independent variables
.fieller.ind <- function(auc1, auc2, var1, var2, df, conf.level=0.90){

  alpha <- 1-conf.level
  # 1-alpha confidence interval
  t <- qt(p=1-alpha/2, df)
  lower <- auc1*auc2 - sqrt((auc1*auc2)^2 - (auc1^2 - t^2*var1)*(auc2^2 - t^2*var2))
  upper <- auc1*auc2 + sqrt((auc1*auc2)^2 - (auc1^2 - t^2*var1)*(auc2^2 - t^2*var2))
  lower <- lower / (auc2^2 - t^2*var2)
  upper <- upper / (auc2^2 - t^2*var2)
 
  # check condition
  res <- NA
  tmax <- (auc1^2*var2 + auc2^2*var1)/(var1*var2)
  if(0 < t^2 & t^2 < tmax^2){res <- c(lower, upper)}
  return(res)
}


# Fieller confidence interval for dependent variables
.fieller.general <- function(auc1, auc2, var1, var2, covar, df, conf.level=0.90){

  alpha <- 1-conf.level
  if(is.na(df)){
    c <- qnorm(1-alpha/2)
  }else{
    c <- qt(1-alpha/2, df=df)
  }

  k <- c^2 * var2/auc2^2
  delta <- auc1/auc2  
  part1 <- delta + (k/(1-k))*(delta - covar/var2)  
  part2 <- c/(auc2*(1-k)) * sqrt(var1 - 2*delta*covar + delta^2*var2 - k * (var1 - (covar^2 / var2)))
  
  upper <- part1 + part2
  lower <- part1 - part2  
  res <- c(lower, upper)
  return(res)
}

# function to modify a dataframe into the list format for batch design
.formattobatch <- function(data){

  data <- data[order(data$time, data$id),]
  J <- length(unique(data$time))
  n <- length(unique(data$id))
  concmat<- matrix(NA,nrow=J,ncol=n)
  rownames(concmat) <- sort(unique(data$time))
  colnames(concmat) <- sort(unique(data$id))
  for(i in 1:nrow(data)){
    concmat[paste(data[i,'time']),paste(data[i,'id'])] <- data[i,'conc']
  }
  #incidence matrix
  D <-matrix(0,nrow=J,ncol=n)
  D[!is.na(concmat)] <- 1

  unique_batches <- t(unique(t(D)))
  # number of batches
  nb <- ncol(unique_batches)
  b <- matrix(NA,ncol=n,nrow=1)
  colnames(b)<- colnames(concmat)
  for( i in 1:nb){
    b[colSums(D==unique_batches[,i])==J] <- i
  }

  for(i in 1:nrow(data)){
    data[i,'B'] <-  b[colnames(b)==data[i,'id']]
  }

  return(split(data,data$B))

}

# Function to convert lists to matrix with one row per subject and time
.batchtoFDAsend<-function(conc,time,group){

  ## sample size within each batch
  ns <- unlist(lapply(time,length))/unlist(lapply(lapply(time,unique),length))

  # create list of ids
  i<-1
  ids <- time

  for (b in 1:length(conc)){
    ids[[b]]<-rep(seq(i,cumsum(ns)[b]),length(unique(time[[b]])))
    i<-i+ns[b]
  }

  if(!is.null(group)){
    data<-data.frame(id=unlist(ids),time=unlist(time),conc=unlist(conc),group=unlist(group))
  }else{
    data<-data.frame(id=unlist(ids),time=unlist(time),conc=unlist(conc))
  }
  data <- data[order(data$time, data$id),]
  return(data)
}

# function to modify a set of lists into a datamatrix and an incidence matrix
.batchtomatrix<-function(conc,time){

  ## sample size within each batch
  ns <- unlist(lapply(time,length))/unlist(lapply(lapply(time,unique),length))

  # create list of ids
  i<-1
  ids <- time

  for (b in 1:length(conc)){
    #sort time and conc so that same times are next to each other
    conc[[b]] <- conc[[b]][order(time[[b]])]
    time[[b]] <- sort(time[[b]])
    ids[[b]]<-rep(seq(i,cumsum(ns)[b]),length(unique(time[[b]])))
    i<-i+ns[b]
  }

  data<-data.frame(id=unlist(ids),time=unlist(time),conc=unlist(conc))
  data <- data[order(data$time, data$id),]
  J <- length(unique(data$time))
  n <- length(unique(data$id))
  concmat<- matrix(NA,nrow=J,ncol=n)

  rownames(concmat) <- sort(unique(data$time))
  colnames(concmat) <- sort(unique(data$id))
  for(i in 1:nrow(data)){
    concmat[paste(data[i,'time']),paste(data[i,'id'])] <- data[i,'conc']
  }

  return(concmat)
}

#obj coming from test()
.resampling.test <- function(obj, theta, nsample=1000, alternative="two.sided"){

## find out what parameters were used.
  if(nrow(obj$est)>1) {
    parm <- 'nca'
    fun <- nca
  }else{
    perm.dist <- rep(NA, nsample) ## storage for permutation distribution
    if(substr(rownames(obj$est)[1],1,5)=='ratio'){
      parm <- 'eqv'
      fun <- eqv
    }else{
      parm <- 'auc'
      fun <- auc
    }
  }

  perm.dist <- matrix(NA, ncol=nsample,nrow=nrow(obj$est)) ## storage for permutation distribution
  stat.obs <- (obj$est[,1]-theta)/obj$CIs[,2] # observed statistic

  if(!is.null(obj$group) && parm=="auc"){ # use permutation test    
    if(obj$design=="ssd"){ # permutation test based on ptest.ssd from Version 1.01

      data <- data.frame(conc=obj$conc, time=obj$time, group=factor(obj$group))
      data <- data[order(data$time, data$group),]

      intern <- data # internal object
      for(i in 1:nsample){
        intern$ran <- runif(nrow(intern))
        intern <- intern[order(intern$time, intern$ran),]
        intern$group <- data$group
        temp <- auc.ssd(data=intern,nsample=0)
        perm.dist[1,i] <- (temp$est[1,1]-theta)/temp$CIs[1,2]
      }
      switch(alternative, 
        "less"={p.value <- length(subset(perm.dist[1,], stat.obs <= perm.dist[1,]))/nsample},	
        "greater"={p.value <- length(subset(perm.dist[1,], stat.obs >= perm.dist[1,]))/nsample},
        "two.sided"={p.value <- length(subset(perm.dist[1,], abs(perm.dist[1,])>abs(stat.obs)))/nsample},
      )
      return(p.value)
    }else{ # permutation test for batch/complete data design
   
      B <- length(obj$conc)
      conc <- obj$conc
      for(i in 1:nsample){
        temp <- NULL       
        for(j in 1:B){
          temp <- data.frame(conc=obj$conc[[j]],time=obj$time[[j]])
          n <- as.numeric(table(temp$time)[1])
          temp <- temp[order(temp$time),]
          temp$ran <- runif(n)
          temp <- temp[order(temp$time, temp$ran),]
          conc[[j]] <- temp$conc
        }
        temp <- auc.batch(conc=conc,time=lapply(obj$time,sort),obj$group,method="z")
        perm.dist[1,i] <- (temp$est[1,1]-theta)/temp$CIs[1,2] 
      }
      switch(alternative, 
        "less"={p.value <- length(subset(perm.dist[1,], stat.obs <= perm.dist[1,]))/nsample},	
        "greater"={p.value <- length(subset(perm.dist[1,], stat.obs >= perm.dist[1,]))/nsample},
        "two.sided"={p.value <- length(subset(perm.dist[1,], abs(perm.dist[1,])>abs(stat.obs)))/nsample},
      )
      return(p.value)      
    }

  }else{ # use one-sample bootstrap test

    boott <- .boot.t.stat(obj, fun=fun, nsample=nsample)

    stat <- (theta - obj$est)/obj$CIs[obj$CIs[,'method']==obj$CIs[1,'method'],2]

    p <- seq(0,0.5,0.0001)
    tstar_alpha <- t(apply(boott,1,quantile,c(p,1-p),type=5))

    switch(alternative, 
      "less"={p.value <- (p[apply(abs(tstar_alpha[,(length(p)+1):(length(p)*2)]-rep(stat,length(p))),1,which.min)])*2},	
      "greater"={p.value <- (p[apply(abs(tstar_alpha[,1:length(p)]-rep(stat,length(p))),1,which.min)])*2},
      "two.sided"={p.value <- (c(p,p)[apply(abs(tstar_alpha-rep(stat,length(p)*2)),1,which.min)])*2},
    )  
    return(p.value)

  }

}

# fun is function call for parameter estimation
.boot.t.stat <- function(obj,fun, nsample){

  stat <- matrix(NA, ncol=nsample,nrow=nrow(obj$est)) ## storage for permutation distribution
  if(obj$design=='ssd'){

    group <- obj$group
    if(is.null(obj$group)) group <- rep(1,length(obj$conc)) 
    n <-   as.numeric(table(obj$time))
    data <- data.frame(conc=obj$conc, time=obj$time, group=factor(group))
    data <- data[order(data$time, data$group),]
    intern <- data
    for(i in 1:nsample){

      inds <- as.vector(unlist(lapply(split(1:nrow(data),obj$time),sample,n,replace=TRUE)))
      intern$conc <- data$conc[inds]
      intern$group <- data$group
      if(is.null(obj$dose)){
        temp <- fun(data=intern,method="z",design='ssd')
      }else{
        temp <- fun(data=intern,method="z",dose=obj$dose, design='ssd')
      }
      stat[,i] <- (temp$est[,1]-obj$est[,1])/temp$CIs[,2] # bootstrap-t-vector
    }
  }else{  # batch design

    B <- length(obj$conc)
    conc <- obj$conc
    for(i in 1:nsample){
      temp <- data <- NULL       
      for(j in 1:B){
        temp <- data.frame(conc=obj$conc[[j]],time=obj$time[[j]])
        n <- as.numeric(table(temp$time))
        inds <- as.vector(unlist(lapply(split(1:nrow(temp),temp$time),sample,n,replace=TRUE)))
        conc[[j]] <- temp$conc[inds]
      }
      temp <- fun(conc=conc,time=obj$time,obj$group,method="z", design="batch")
      stat[,i] <- (temp$est[,1]-obj$est[,1])/temp$CIs[,2] # bootstrap-t-vector
    }
  }

  return(stat)

}


