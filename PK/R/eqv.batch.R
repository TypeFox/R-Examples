eqv.batch <- function(conc, time, group, dependent=FALSE, method=c("fieller", "z", "boott"), conf.level=0.90, nsample=1000, data){

  ### function to find bootstrap-t interval
  boott <- function(conc1, conc2, time1, time2, asym.var, dependent, conf.level, nsample){	
    # get estimates
    obsv.auc1 <- auc.batch(conc1, time1, method='z')$CIs[1,1:2]
    obsv.auc2 <- auc.batch(conc2, time2, method='z')$CIs[1,1:2]
    obsv.est <- obsv.auc1$est / obsv.auc2$est
    obsv.var <- asym.var

    n <- as.numeric(table(time1[[1]])[1])			
    m <- as.numeric(table(time2[[1]])[1])
    B <- length(time1)
    boot.conc1<-vector('list',B)
    boot.conc2<-vector('list',B)	
    # bootstrap distribution
    boot.stat <- rep(NA,nsample)
    for(j in 1:nsample){
      for (i in 1:B){
        boot.conc1[[i]]<-as.vector(matrix(conc1[[i]],nrow=n)[sample(1:n,replace=TRUE),])
        boot.conc2[[i]]<-as.vector(matrix(conc2[[i]],nrow=m)[sample(1:m,replace=TRUE),])
      }
      boot.auc1 <- auc.batch(boot.conc1, time1, method="t",  alternative="two.sided", conf.level=conf.level, nsample=0)$CIs[1,1:2]
      boot.auc2 <- auc.batch(boot.conc2, time2, method="t",  alternative="two.sided", conf.level=conf.level, nsample=0)$CIs[1,1:2]
      v11 <- boot.auc1$stderr^2
      v22 <- boot.auc2$stderr^2
      boot.est <- boot.auc1$est / boot.auc2$est
      if(dependent){
        boot.dconc <- NULL
        for(i in 1:length(time1)){
          boot.dconc[[i]] <- boot.conc1[[i]]-boot.conc2[[i]]
        }
        boot.dauc <- auc.batch(conc=boot.dconc, time=time1, group=NULL, method="t",  alternative="two.sided", conf.level=conf.level, nsample=0)
        v12 <- 0.5*(v11+v22-boot.dauc$CIs[1,2]^2)
        boot.var <- v11/a2^2+v22*a1^2/a2^4-2*v12*a1/a2^3
      }else{
        boot.var <- v11/a2^2+a1^2/a2^4*v22
      }
      boot.stat[j] <- (boot.est - obsv.est) / sqrt(boot.var)
    }

    alpha <- 1-conf.level
    t.lb <- quantile(boot.stat, probs=c(alpha/2),   method=5, na.rm=TRUE)
    t.ub <- quantile(boot.stat, probs=c(1-alpha/2), method=5, na.rm=TRUE)
    base <- data.frame(est=obsv.est, t.lb=t.lb, t.ub=t.ub)	
    base$lower <- base$est - base$t.ub*sqrt(obsv.var)
    base$upper <- base$est - base$t.lb*sqrt(obsv.var)
    return(c(base$lower, base$upper))
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
  method <- sort(method)
  if(!is.null(group)){
    if(length(unlist(conc))!=length(unlist(group))){stop('different length of concentration and grouping lists')}
    if(nlevels(as.factor(unlist(group))) > 2){stop("limited for comparison of 2 groups")}
    if(nlevels(as.factor(unlist(group))) == 1) {group <- NULL}
  }

  # handle input parameters
  if(any(is.na(unlist(conc)))){stop('Missing values are not permitted in batch designs')}
  if (!is.list(time) || !is.list(conc)) stop('Both time and concentration need to be a list')
  if (!is.list(group)) stop('The parameter group needs to be a list')
  if (length(time) != length(conc)) stop('Number of batches defined in time and conc are different')
  if (length(time) != length(group)) stop('Number of batches defined in group does not match time and conc')

  ##CHECK FOR SAMPLE SIZE TO BE EQUAL

  conc1 <- NULL
  conc2 <- NULL
  time1 <- NULL
  time2 <- NULL
  grps<-unique(unlist(group))
  if(length(grps)!=2)stop('Number of groups is unequal 2')
  B<-length(conc)
  for(i in 1:B){
    conc1[[i]] <- conc[[i]][group[[i]]==grps[1]]
    conc2[[i]] <- conc[[i]][group[[i]]==grps[2]]
    time1[[i]] <- time[[i]][group[[i]]==grps[1]]
    time2[[i]] <- time[[i]][group[[i]]==grps[2]]
  }
  names(conc1)<-names(conc)
  names(conc2)<-names(conc)
  names(time1)<-names(time)
  names(time2)<-names(time)

  auc1 <- auc.batch(conc1, time1, method="t",  alternative="two.sided", conf.level=conf.level, nsample=0)
  auc2 <- auc.batch(conc2, time2, method="t",  alternative="two.sided", conf.level=conf.level, nsample=0)
 
  a1 <- auc1$CIs[1,1]
  a2 <- auc2$CIs[1,1]
  est <- a1/a2
  v11 <- auc1$CIs[1,2]^2
  v22 <- auc2$CIs[1,2]^2

  if(dependent){
    #compute difference of AUCs 
    dconc <- NULL
    for(i in 1:length(time)){
      dconc[[i]] <- conc1[[i]]-conc2[[i]]
    }
    names(dconc) <- names(conc)
    dauc <- auc.batch(conc=dconc, time=time1, group=NULL, method="t",  alternative="two.sided", conf.level=conf.level, nsample=0)
    v12 <- 0.5*(v11+v22-dauc$CIs[1,2]^2)
    asym.var <- v11/a2^2+v22*a1^2/a2^4-2*v12*a1/a2^3
  }else{
    asym.var <- v11/a2^2+a1^2/a2^4*v22
  }

  res <- NULL
  df <- NULL
  if(any(method=='boott')){
    if(nsample<2) stop('nsample needs to be >1 for method boott')
    res <- rbind(res,boott(conc1=conc1, conc2=conc2, time1=time1, time2=time2, asym.var=asym.var, dependent=dependent, conf.level=conf.level, nsample=nsample))
    df <- rbind(df,NA)
  }
  if(any(method=='fieller')){
    if(dependent){
      res <- rbind(res,.fieller.general(auc1=a1 , auc2=a2, var1=v11, var2=v22, covar=v12, df=dauc$CIs[1,5],conf.level=conf.level))
      df <- rbind(df,dauc$CIs[1,5])
    }else{
      n <- as.numeric(table(time1[[1]])[1])
      m <- as.numeric(table(time2[[1]])[1])
      df.ind<-(v11+(a1/a2)^2*v22)^2/(v11^2/(B*(n-1))+(a1/a2)^4*v22^2/(B*(m-1))) 
      res <- rbind(res,.fieller.ind(auc1=a1, auc2=a2, var1=v11, var2=v22, df=df.ind, conf.level=conf.level))
      df <- rbind(df,df.ind)
    }
  }
  if(any(method=='z')){
    alpha <- 1-conf.level
    res <- rbind(res,a1/a2+qnorm(c(alpha/2,1-alpha/2))*sqrt(asym.var))
    df <- rbind(df,NA)
  }
  rownames(res) <- method
  rownames(df) <- method
  colnames(res) <- c('lower','upper')

  out <- NULL
  out$est <- matrix(as.double(est),ncol=1)
  if(dependent){
    rownames(out$est) <- 'ratio of dependent AUCs to tlast'
  }else{
    rownames(out$est) <- 'ratio of independent AUCs to tlast'
  }
  colnames(out$est) <- 'est'

  out$CIs<-data.frame(est=rep(est,length(method)), stderr=rep(sqrt(asym.var),length(method)), lower=res[,1], upper=res[,2], df=df ,method=method)
  if(dependent){
    rownames(out$CIs) <- paste(conf.level*100,'% CI using a ', method,'-interval for the ratio of dependent AUCs to tlast', sep='')
  }else{
    rownames(out$CIs) <- paste(conf.level*100,'% CI using a ', method,'-interval for the ratio of independent AUCs to tlast', sep='')
  }
  out$design<-"batch"
  out$conf.level <- conf.level
  out$conc <- conc
  out$time <- time
  out$group <- group
  class(out)<-"PK"
  return(out)
}
