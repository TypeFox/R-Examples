auc.batch <- function(conc, time, group=NULL, method=c("t", "z", "boott"),  alternative=c("two.sided", "less", "greater"), conf.level=0.95, nsample=1000, data){
 
  # function for bootstrap resampling
  resample <- function(w, conc, time, group, alpha, nsample){
    B<-length(time)
    obsv.parm <- parms(w=w, conc=conc, time=time, group=group)

    boot.conc<-vector('list',B)
    boot.group<-vector('list',B)
    boot.time <-vector('list',B)
    stat <- rep(NA, nsample) 
    for(j in 1:nsample){
       ## generate bootstrap data
       if(!is.null(group)){
         grps<-unique(unlist(group))
         n<-length(conc[[1]][group[[1]]==grps[1]])/length(unique(time[[1]]))
         m<-length(conc[[1]][group[[1]]==grps[2]])/length(unique(time[[1]]))
         for (i in 1:B){
           boot.conc[[i]]<-c(as.vector(matrix(conc[[i]][group[[i]]==grps[1]],nrow=n)[sample(1:n,replace=TRUE),]),as.vector(matrix(conc[[i]][group[[i]]==grps[2]],nrow=m)[sample(1:m,replace=TRUE),]))
           boot.group[[i]]<-c(rep(1,table(group[[i]])[1]),rep(2,table(group[[i]])[2]))
           boot.time[[i]] <- c(rep(unique(time[[i]]),each=n),rep(unique(time[[i]]),each=m))
         }
       }else{
         n<-length(conc[[1]])/length(unique(time[[1]]))
         for (i in 1:B){
           boot.conc[[i]]<-as.vector(matrix(conc[[i]],nrow=n)[sample(1:n,replace=TRUE),])
         }
         boot.time <- time
         boot.group<-NULL
       }

       boot.parm <- parms(w=w, conc=boot.conc, time=boot.time, group=boot.group)
       stat[j] <- (boot.parm$est- obsv.parm$est) / sqrt(boot.parm$var)
    }
    lower <- quantile(stat, alpha/2, method=5)
    upper <- quantile(stat, 1-alpha/2, method=5)
    return(c(upper, lower*(-1)))
  }

  # function to calculate AUC and V(AUC)
  parms <- function(w, conc, time, group){

    single.AUC <- function(w, conc, time){

      concmat <- .batchtomatrix(conc,time)
      J <- nrow(concmat)
      ### total number of subjects
      n.tot <- ncol(concmat)
      #incidence matrix
      Inc <-matrix(0,nrow=J,ncol=n.tot)
      Inc[!is.na(concmat)] <- 1
      
      delta <- matrix(NA,ncol=n.tot,nrow=J)
      for(i in 1:n.tot){
        delta[,i] <- sum(colSums(Inc[,i]==Inc)==J)/rowSums(Inc)
      }
      delta[Inc==0]<- NA

  ## compute AUC
  # find which batch each observation belongs to
      unique_batches <- t(unique(t(Inc)))
      nb <- ncol(unique_batches)
      b <- rep(NA,n.tot)
      for( i in 1:nb){
        b[colSums(Inc==unique_batches[,i])==J] <- i
      }

      delta[is.na(delta)] <- 0
      AUC <- sum(tapply(colSums(delta*matrix(1,ncol=n.tot,nrow=J)*w*concmat,na.rm=TRUE),b,mean))  # partial AUC for each individual

      ### number of observation following each batch
      n <- table(b)
      R <- Inc %*% t(Inc)
      rvec <- rowSums(Inc)

      ## average per timepoint
      y.mean <- rowMeans(concmat, na.rm=TRUE)
      ## deviation from mean
      y.dif <- apply(concmat,2,"-",y.mean)

      ## covariance 
      cov.jh <- matrix(NA,nrow(concmat),nrow(concmat))
      diag(cov.jh)<- diag(var(t(concmat),na.rm=TRUE,use="pairwise.complete.obs"))
      for (j1 in 1:(nrow(concmat))) {
        for(j2 in (j1):nrow(concmat)) {
          r <- R[j1,j2]
          cov.jh[j1,j2] <- sum(y.dif[j1,]*y.dif[j2,],na.rm=TRUE)/((r-1)+(1-r/rvec[j1])*(1-r/rvec[j2]))
          if(!is.nan(cov.jh[j1,j2]) & abs(cov.jh[j1,j2]) > sqrt(cov.jh[j1,j1]*cov.jh[j2,j2])){
            cov.jh[j1,j2] <- sign(cov.jh[j1,j2])*sqrt(cov.jh[j1,j1]*cov.jh[j2,j2])
          }
        }
      }

      ind <-lower.tri(cov.jh)
      cov.jh[ind]<-t(cov.jh)[ind]

      cov.jh[is.nan(cov.jh)] <- 0
      vAUC <- 0
      temp <- rep(NA,n.tot)
      for(i in 1:n.tot){
        help <- 0
        for(j1 in 1:J){
          help <- help + sum(delta[j1,i]*delta[,i]*w[j1]*w*cov.jh[j1,])
        }
        temp[i] <- 1/n[b[i]]^2*help                
      }
      # vector of variance of partial AUCs
      vpAUC <- tapply(temp,b,sum)
      vAUC <- sum(vpAUC)
      df <- vAUC^2/sum(vpAUC^2/((n-1)))

      res <- list(var=vAUC, est=AUC, df=df)
      return(res)

    }

    B<-length(time) ### number of batches

    ## difference of two AUC's
    if(!is.null(group)){
 
      conc1 <- conc2 <- time1<-time2 <- vector('list',B)
      grps<-unique(unlist(group))
      
      # seperating groups
      for(i in 1:B){
        conc1[[i]]<-conc[[i]][group[[i]]==grps[1]]
        time1[[i]]<-time[[i]][group[[i]]==grps[1]]
        conc2[[i]]<-conc[[i]][group[[i]]==grps[2]]
        time2[[i]]<-time[[i]][group[[i]]==grps[2]]
      }
      
      res1 <- single.AUC(w,conc1,time1)
      res2 <- single.AUC(w,conc2,time2)

      res <- list(var=res1$var+res2$var, est=res1$est-res2$est, df=res1$df+res2$df)
      return(res)
    }else{

      ## single AUC      
    
      return(single.AUC(w,conc,time))

    }
  }

  # check input parameters  
  if(!missing(data)){
    cnames <- colnames(data)
    if(!any(cnames=='id')){stop("data does not contain a variable id")}
    if(!any(cnames=='conc')){stop("data does not contain a variable conc")}
    if(!any(cnames=='time')){stop("data does not contain a variable time")}
    temp <- .formattobatch(data)
    conc <- time <- group <- NULL
    for(i in 1:length(temp)){
      conc[[i]] <- temp[[i]]$conc 
      time[[i]] <- temp[[i]]$time
      if(any(names(temp[[1]])=='group')){
        group[[i]] <- temp[[i]]$group
      }
    }
  }  

  method <- match.arg(method,several.ok=TRUE)
  alternative <- match.arg(alternative)
  if(!is.null(group)){
    if(length(unlist(conc))!=length(unlist(group))){stop('different length of concentration and grouping lists')}
    if(nlevels(as.factor(unlist(group))) > 2){stop("limited for comparison of 2 groups")}
    if(nlevels(as.factor(unlist(group))) == 1) {group <- NULL}
  }
  
  # handle input parameters
  if(any(is.na(unlist(conc)))){stop('Missing values are not permitted in batch designs')}
  if (!is.list(time) || !is.list(conc)) stop('Both time and concentration need to be a list')
  if (!is.null(group) && !is.list(group)) stop('The parameter group needs to be a list or NULL')
  if (length(time) != length(conc)) stop('Number of batches defined in time and conc are different')
  if (length(time) != length(group) && !is.null(group)) stop('Number of batches defined in group does not match time and conc')


  n<-as.numeric(table(time[[1]])[1])
#  if(sum(as.vector(table(unlist(time)))!=rep(n,times=length(unique(unlist(time)))))>0) stop('Number of observations differs between batches')
#  if (sum(table(unlist(lapply(time,unique)))>1)>0) stop('Time points in each batch are not unique')

  grpfact <- levels(as.factor(unlist(group)))
  # specify alpha error
  alpha <- 1-conf.level

  if(is.null(group)){
    grp <- lapply(lapply(conc,'>',-Inf),as.numeric)
  }else{
    grp <- group
  }

  # sort concentrations by time order
  for(i in 1:length(time)){
    data <- data.frame(conc=conc[[i]], time=time[[i]], group=grp[[i]])
    data <- data[order(data$time), ]
    conc[[i]]<-data$conc
    time[[i]]<-data$time
    grp[[i]] <- data$group
  }

  if(!is.null(group)){
    group<-grp
  }  

  # eliminates multiples in time points
  time.unique <- lapply(time,unique)
  # calculate weights and observed parameters
  w <- .weight(time.unique)  

  # calculate confidence intervals

  obsv.parm <- parms(w=w, conc=conc, time=time, group=group)
  if(alternative %in% c('less', 'greater')){alpha <- alpha*2}    
  z <- NULL
  
  method <- sort(method)  
  if (any(method=="boott")){
    z <- rbind(z,resample(w=w, conc=conc, time=time, group=group, alpha=alpha, nsample=nsample))
  }

  if (any(method=="t")){
    z <- rbind(z,rep(qt(1-alpha/2, df=obsv.parm$df),2))
  }

  if (any(method=="z")) {
    z <- rbind(z,rep(qnorm(1-alpha/2),2))
  }
       
  est <- sum(obsv.parm$est[1], -obsv.parm$est[2], na.rm=TRUE)    
  lower <- est - sqrt(sum(obsv.parm$var))*z[,1]
  upper <- est + sqrt(sum(obsv.parm$var))*z[,2]
    
  switch(alternative,
    "less"={upper <- Inf},
    "greater"={lower <- -Inf},
    "two.sided"={},
  )

  df <- rep(NA, as.double(length(lower)))
  if(any(method=="t")){
    if(any(method=='z')){
      df=c(rep(NA, as.double(length(lower)-2)), obsv.parm$df, NA)
    }else{
      df=c(rep(NA, as.double(length(lower)-1)), obsv.parm$df)
    }
  } 
  res <- NULL
  res$est <- matrix(as.double(est),ncol=1)
  if(is.null(group)){
    rownames(res$est) <- 'AUC to tlast'
  }else{
    rownames(res$est) <- 'difference of AUCs to tlast'
  }
  colnames(res$est) <- 'est'
  res$design<-"batch"
  res$CIs<-data.frame(est=est, stderr=sqrt(sum(obsv.parm$var)), lower=lower, upper=upper, df=df,method=method)
  rownames(res$CIs) <- paste(conf.level*100,'% CI using a ', method,' distribution for AUC to tlast', sep='')
  res$conf.level <- conf.level
  res$conc <- conc
  res$time <- time
  res$group <- group
  class(res)<-"PK"
  return(res)
}


