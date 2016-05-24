nma.ab.bin <-
function(s.id,t.id,event.n,total.n,data,trtname,param=c("AR","LOR","LRR","RD","rank.prob"),model="het_cor",prior.type,a=0.001,b=0.001,c=10,higher.better=FALSE,digits=4,n.adapt=5000,n.iter=100000,n.burnin=floor(n.iter/2),n.chains=3,n.thin=max(1,floor((n.iter-n.burnin)/100000)),conv.diag=FALSE,trace=NULL,dic=FALSE,postdens=FALSE,mcmc.samples=FALSE){
  ## check the input parameters
  options(warn=1)
  if(missing(s.id)) stop("need to specify study id.")
  if(missing(t.id)) stop("need to specify treatment.")
  if(missing(event.n)) stop("need to specify event number.")
  if(missing(total.n)) stop("need to specify total number.")
  if(!missing(data)){
    s.id<-eval(substitute(s.id),data,parent.frame())
    t.id<-eval(substitute(t.id),data,parent.frame())
    event.n<-eval(substitute(event.n),data,parent.frame())
    total.n<-eval(substitute(total.n),data,parent.frame())
  }
  if(length(s.id)!=length(t.id) | length(t.id)!=length(event.n) | length(event.n)!=length(total.n) | length(total.n)!=length(s.id)){
    stop("s.id, t.id, event.n, and total.n have different lengths.")
  }
  if(!all(event.n<=total.n)) stop("total number must be greater than event number.")
  if(!all(total.n>0)) stop("total number must be positive.")
  if(!all(event.n>=0)) stop("event number must be non-negative.")
  if(!all(event.n%%1==0) | !all(total.n%%1==0)) warning("at least one event number or total number is not integer.")
  if(!is.element(model,c("hom_eqcor","het_eqcor","het_cor"))) stop("model should be specified as \"hom_eqcor\", \"het_eqcor\", or \"het_cor\".")

  if(any(is.na(s.id))|any(is.na(t.id))|any(is.na(event.n))|any(is.na(total.n))){
    dat<-cbind(s.id,t.id,event.n,total.n)
    s.id<-s.id[complete.cases(dat)]
    t.id<-t.id[complete.cases(dat)]
    event.n<-event.n[complete.cases(dat)]
    total.n<-total.n[complete.cases(dat)]
    cat("NA is not allowed in the input data set;\n")
    cat("the rows containing NA are removed.\n")
  }

  ## make ids continuous
  s.id.o<-s.id
  t.id.o<-t.id
  s.label<-sort(unique(s.id.o))
  t.label<-sort(unique(t.id.o))
  nstudy<-length(s.label) # total number of studies
  ntrt<-length(t.label) # total number of treatments
  len<-length(s.id)
  s.id<-numeric(nstudy)
  for(i in 1:nstudy) {s.id[which(s.id.o==s.label[i])]<-i}
  t.id<-numeric(ntrt)
  for(i in 1:ntrt) {t.id[which(t.id.o==t.label[i])]<-i}

  if(missing(trtname)){
    if(is.numeric(t.id.o)) {trtname<-paste("Trt",t.label,sep="")}
    if(is.character(t.id.o)) {trtname<-t.label}
  }
  if(length(trtname)!=length(unique(t.id))) stop("the number of treatment names does not match for specified treatment id.")
  if(missing(prior.type)) prior.type<-ifelse(model=="het_cor","invwishart","unif")

  ## jags model
  if(model=="hom_eqcor"){
    modelstring<-model.binary.hom.eqcor(prior.type,is.element("rank.prob",param))
  }
  if(model=="het_eqcor"){
    modelstring<-model.binary.het.eqcor(prior.type,is.element("rank.prob",param))
  }
  if(model=="het_cor"){
    modelstring<-model.binary.het.cor(prior.type,is.element("rank.prob",param))
  }

  ## jags data
  if(model == "hom_eqcor"| model == "het_eqcor"){
    if(prior.type == "unif"){
      if(is.element("rank.prob",param)) data.jags<-list(s=s.id,t=t.id,r=event.n,totaln=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),c=c,higher.better=higher.better)
      if(!is.element("rank.prob",param)) data.jags<-list(s=s.id,t=t.id,r=event.n,totaln=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),c=c)
    }
    if(prior.type == "invgamma"){
      if(is.element("rank.prob",param)) data.jags<-list(s=s.id,t=t.id,r=event.n,totaln=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),a=a,b=b,higher.better=higher.better)
      if(!is.element("rank.prob",param)) data.jags<-list(s=s.id,t=t.id,r=event.n,totaln=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),a=a,b=b)
    }
  }
  if(model=="het_cor"){
    if(prior.type == "invwishart"){
      I <- diag(ntrt)
      if(is.element("rank.prob",param)) data.jags<-list(s=s.id,t=t.id,r=event.n,totaln=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),I=I,higher.better=higher.better)
      if(!is.element("rank.prob",param)) data.jags<-list(s=s.id,t=t.id,r=event.n,totaln=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),I=I)
    }
    if(prior.type == "chol"){
      if(is.element("rank.prob",param)) data.jags<-list(s=s.id,t=t.id,r=event.n,totaln=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),c=c,higher.better=higher.better)
      if(!is.element("rank.prob",param)) data.jags<-list(s=s.id,t=t.id,r=event.n,totaln=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),c=c)
    }
  }

  ## jags initial value
  rng.seeds<-sample(1000000,n.chains)
  mu.init<-numeric(ntrt)
  for(i in 1:ntrt){
    mu.init[i]<-sum(event.n[t.id==t.id[i]])/sum(total.n[t.id==t.id[i]])
  }
  init.jags<-list(NULL)
  if(model=="hom_eqcor"){
    if(prior.type=="unif"){
      for(ii in 1:n.chains){
        init.jags[[ii]]<-list(mu=qnorm(mu.init),vi=matrix(0,nstudy,ntrt),sigma=c/2,rho=0.5,.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
      }
    }
    if(prior.type=="invgamma"){
      for(ii in 1:n.chains){
        init.jags[[ii]]<-list(mu=qnorm(mu.init),vi=matrix(0,nstudy,ntrt),inv.sig.sq=a/b,rho=0.5,.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
      }
    }
  }
  if(model=="het_eqcor"){
    if(prior.type=="unif"){
      for(ii in 1:n.chains){
        init.jags[[ii]]<-list(mu=qnorm(mu.init),vi=matrix(0,nstudy,ntrt),sigma=rep(c/2,ntrt),rho=0.5,.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
      }
    }
    if(prior.type=="invgamma"){
      for(ii in 1:n.chains){
        init.jags[[ii]]<-list(mu=qnorm(mu.init),vi=matrix(0,nstudy,ntrt),inv.sig.sq=rep(a/b,ntrt),rho=0.5,.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
      }
    }
  }
  if(model=="het_cor"){
    if(prior.type=="invwishart"){
      for(ii in 1:n.chains){
        init.jags[[ii]]<-list(mu=qnorm(mu.init),vi=matrix(0,nstudy,ntrt),T=(ntrt+1)*I,.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
      }
    }
    if(prior.type=="chol"){
      for(ii in 1:n.chains){
        init.jags[[ii]]<-list(mu=qnorm(mu.init),vi=matrix(0,nstudy,ntrt),sigma=rep(c/2,ntrt),psi=matrix(3.1415926/2,ntrt-1,ntrt-1),.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
      }
    }
  }

  ## parameters to be monitored in jags
  if(!is.element("AR",param)) param<-c("AR",param)
  if(!is.null(trace)){
    if(!any(is.element(trace, param))) stop("at least one effect size in argument trace is not specified in argument param.")
  }
  monitor<-param[!is.element(param,c("OR","RR","RD","LOR","LRR"))]
  if(is.element("RR",param)){
    for(ii in 1:ntrt){
      for(jj in 1:ntrt){
        if(ii != jj) monitor<-c(monitor,paste("RR[",ii,",",jj,"]",sep=""))
      }
    }
  }
  if(is.element("RD",param)){
    for(ii in 1:ntrt){
      for(jj in 1:ntrt){
        if(ii < jj) monitor<-c(monitor,paste("RD[",ii,",",jj,"]",sep=""))
      }
    }
  }
  if(is.element("OR",param)){
    for(ii in 1:ntrt){
      for(jj in 1:ntrt){
        if(ii != jj) monitor<-c(monitor,paste("OR[",ii,",",jj,"]",sep=""))
      }
    }
  }
  if(is.element("LRR",param)){
    for(ii in 1:ntrt){
      for(jj in 1:ntrt){
        if(ii < jj) monitor<-c(monitor,paste("LRR[",ii,",",jj,"]",sep=""))
      }
    }
  }
  if(is.element("LOR",param)){
    for(ii in 1:ntrt){
      for(jj in 1:ntrt){
        if(ii < jj) monitor<-c(monitor,paste("LOR[",ii,",",jj,"]",sep=""))
      }
    }
  }

  ## run jags
  cat("Start running MCMC...\n")
  jags.m<-tryCatch.W.E(jags.model(file=textConnection(modelstring),data=data.jags,inits=init.jags,n.chains=n.chains,n.adapt=n.adapt))
  warn.adapt<-jags.m$warning
  jags.m<-jags.m$value
  if(is(warn.adapt,"warning")) cat("Adaptation incomplete; users may increase n.adapt.\n")
  update(jags.m,n.iter=n.burnin)
  jags.out<-coda.samples(model=jags.m,variable.names=monitor,n.iter=n.iter,thin=n.thin)
  smry<-summary(jags.out)
  smry<-cbind(smry$statistics[,c("Mean","SD")],smry$quantiles[,c("2.5%","50%","97.5%")])
  smry<-signif(smry,digits=digits)

  out<-NULL
  out$model<-"Binomial likelihood with probit link."
  AR.id<-grep("AR",rownames(smry))
  AR.stat<-array(paste(format(round(smry[AR.id,"Mean"],digits=digits),nsmall=digits)," (",format(round(smry[AR.id,"SD"],digits=digits),nsmall=digits),")",sep=""),dim=c(ntrt,1))
  colnames(AR.stat)<-"Mean (SD)"
  rownames(AR.stat)<-trtname
  AR.quan<-array(paste(format(round(smry[AR.id,"50%"],digits=digits),nsmall=digits)," (",format(round(smry[AR.id,"2.5%"],digits=digits),nsmall=digits),", ",format(round(smry[AR.id,"97.5%"],digits=digits),nsmall=digits),")",sep=""),dim=c(ntrt,1))
  colnames(AR.quan)<-"Median (95% CI)"
  rownames(AR.quan)<-trtname
  out$AbsoluteRisk<-list(Mean_SD=noquote(AR.stat),Median_CI=noquote(AR.quan))

  if(is.element("OR",param)){
    OR.stat<-OR.quan<-array(NA,dim=c(ntrt,ntrt))
    colnames(OR.stat)<-colnames(OR.quan)<-rownames(OR.stat)<-rownames(OR.quan)<-trtname
    for(i in 1:ntrt){
      OR.stat[i,i]<-OR.quan[i,i]<-"--"
      for(j in 1:ntrt){
        if(i != j){
          OR.ij<-paste("OR[",i,",",j,"]",sep="")
          OR.stat[i,j]<-paste(format(round(smry[OR.ij,"Mean"],digits=digits),nsmall=digits)," (",format(round(smry[OR.ij,"SD"],digits=digits),nsmall=digits),")",sep="")
          OR.quan[i,j]<-paste(format(round(smry[OR.ij,"50%"],digits=digits),nsmall=digits)," (",format(round(smry[OR.ij,"2.5%"],digits=digits),nsmall=digits),", ",format(round(smry[OR.ij,"97.5%"],digits=digits),nsmall=digits),")",sep="")
        }
      }
    }
    out$OddsRatio<-list(Mean_SD=noquote(OR.stat),Median_CI=noquote(OR.quan))
  }

  if(is.element("LOR",param)){
    LOR.stat<-LOR.quan<-array(NA,dim=c(ntrt,ntrt))
    colnames(LOR.stat)<-colnames(LOR.quan)<-rownames(LOR.stat)<-rownames(LOR.quan)<-trtname
    for(i in 1:ntrt){
      LOR.stat[i,i]<-LOR.quan[i,i]<-"--"
      for(j in 1:ntrt){
        if(i < j){
          LOR.ij<-paste("LOR[",i,",",j,"]",sep="")
          LOR.stat[i,j]<-paste(format(round(smry[LOR.ij,"Mean"],digits=digits),nsmall=digits)," (",format(round(smry[LOR.ij,"SD"],digits=digits),nsmall=digits),")",sep="")
          LOR.stat[j,i]<-paste(format(round(-smry[LOR.ij,"Mean"],digits=digits),nsmall=digits)," (",format(round(smry[LOR.ij,"SD"],digits=digits),nsmall=digits),")",sep="")
          LOR.quan[i,j]<-paste(format(round(smry[LOR.ij,"50%"],digits=digits),nsmall=digits)," (",format(round(smry[LOR.ij,"2.5%"],digits=digits),nsmall=digits),", ",format(round(smry[LOR.ij,"97.5%"],digits=digits),nsmall=digits),")",sep="")
          LOR.quan[j,i]<-paste(format(round(-smry[LOR.ij,"50%"],digits=digits),nsmall=digits)," (",format(round(-smry[LOR.ij,"97.5%"],digits=digits),nsmall=digits),", ",format(round(-smry[LOR.ij,"2.5%"],digits=digits),nsmall=digits),")",sep="")
        }
      }
    }
    out$LogOddsRatio<-list(Mean_SD=noquote(LOR.stat),Median_CI=noquote(LOR.quan))
  }

  if(is.element("RR",param)){
    RR.stat<-RR.quan<-array(NA,dim=c(ntrt,ntrt))
    colnames(RR.stat)<-colnames(RR.quan)<-rownames(RR.stat)<-rownames(RR.quan)<-trtname
    for(i in 1:ntrt){
      RR.stat[i,i]<-RR.quan[i,i]<-"--"
      for(j in 1:ntrt){
        if(i != j){
          RR.ij<-paste("RR[",i,",",j,"]",sep="")
          RR.stat[i,j]<-paste(format(round(smry[RR.ij,"Mean"],digits=digits),nsmall=digits)," (",format(round(smry[RR.ij,"SD"],digits=digits),nsmall=digits),")",sep="")
          RR.quan[i,j]<-paste(format(round(smry[RR.ij,"50%"],digits=digits),nsmall=digits)," (",format(round(smry[RR.ij,"2.5%"],digits=digits),nsmall=digits),", ",format(round(smry[RR.ij,"97.5%"],digits=digits),nsmall=digits),")",sep="")
        }
      }
    }
    out$RelativeRisk<-list(Mean_SD=noquote(RR.stat),Median_CI=noquote(RR.quan))
  }

  if(is.element("LRR",param)){
    LRR.stat<-LRR.quan<-array(NA,dim=c(ntrt,ntrt))
    colnames(LRR.stat)<-colnames(LRR.quan)<-rownames(LRR.stat)<-rownames(LRR.quan)<-trtname
    for(i in 1:ntrt){
      LRR.stat[i,i]<-LRR.quan[i,i]<-"--"
      for(j in 1:ntrt){
        if(i < j){
          LRR.ij<-paste("LRR[",i,",",j,"]",sep="")
          LRR.stat[i,j]<-paste(format(round(smry[LRR.ij,"Mean"],digits=digits),nsmall=digits)," (",format(round(smry[LRR.ij,"SD"],digits=digits),nsmall=digits),")",sep="")
          LRR.stat[j,i]<-paste(format(round(-smry[LRR.ij,"Mean"],digits=digits),nsmall=digits)," (",format(round(smry[LRR.ij,"SD"],digits=digits),nsmall=digits),")",sep="")
          LRR.quan[i,j]<-paste(format(round(smry[LRR.ij,"50%"],digits=digits),nsmall=digits)," (",format(round(smry[LRR.ij,"2.5%"],digits=digits),nsmall=digits),", ",format(round(smry[LRR.ij,"97.5%"],digits=digits),nsmall=digits),")",sep="")
          LRR.quan[j,i]<-paste(format(round(-smry[LRR.ij,"50%"],digits=digits),nsmall=digits)," (",format(round(-smry[LRR.ij,"97.5%"],digits=digits),nsmall=digits),", ",format(round(-smry[LRR.ij,"2.5%"],digits=digits),nsmall=digits),")",sep="")
        }
      }
    }
    out$LogRelativeRisk<-list(Mean_SD=noquote(LRR.stat),Median_CI=noquote(LRR.quan))
  }

  if(is.element("RD",param)){
    RD.stat<-RD.quan<-array(NA,dim=c(ntrt,ntrt))
    colnames(RD.stat)<-colnames(RD.quan)<-rownames(RD.stat)<-rownames(RD.quan)<-trtname
    for(i in 1:ntrt){
      RD.stat[i,i]<-RD.quan[i,i]<-"--"
      for(j in 1:ntrt){
        if(i < j){
          RD.ij<-paste("RD[",i,",",j,"]",sep="")
          RD.stat[i,j]<-paste(format(round(smry[RD.ij,"Mean"],digits=digits),nsmall=digits)," (",format(round(smry[RD.ij,"SD"],digits=digits),nsmall=digits),")",sep="")
          RD.stat[j,i]<-paste(format(round(-smry[RD.ij,"Mean"],digits=digits),nsmall=digits)," (",format(round(smry[RD.ij,"SD"],digits=digits),nsmall=digits),")",sep="")
          RD.quan[i,j]<-paste(format(round(smry[RD.ij,"50%"],digits=digits),nsmall=digits)," (",format(round(smry[RD.ij,"2.5%"],digits=digits),nsmall=digits),", ",format(round(smry[RD.ij,"97.5%"],digits=digits),nsmall=digits),")",sep="")
          RD.quan[j,i]<-paste(format(round(-smry[RD.ij,"50%"],digits=digits),nsmall=digits)," (",format(round(-smry[RD.ij,"97.5%"],digits=digits),nsmall=digits),", ",format(round(-smry[RD.ij,"2.5%"],digits=digits),nsmall=digits),")",sep="")
        }
      }
    }
    out$RiskDifference<-list(Mean_SD=noquote(RD.stat),Median_CI=noquote(RD.quan))
  }

  if(is.element("rank.prob",param)){
    rank.prob.id<-grep("rank.prob",rownames(smry))
    rank.prob.stat<-array(format(round(smry[rank.prob.id,"Mean"],digits=4),nsmall=4),dim=c(ntrt,ntrt))
    colnames(rank.prob.stat)<-paste("rank",1:ntrt,sep="")
    rownames(rank.prob.stat)<-trtname
    out$TrtRankProb<-noquote(rank.prob.stat)
  }

  if(conv.diag){
    cat("Start calculating MCMC convergence diagnostic statistics...\n")
    conv.out<-gelman.diag(jags.out,multivariate=FALSE)
    conv.out<-conv.out$psrf
    if(is.element("rank.prob",param)){
      rank.prob.id<-grep("rank.prob",rownames(conv.out))
      conv.out<-conv.out[-rank.prob.id,]
    }
    write.table(conv.out,"ConvergenceDiagnostic.txt",row.names=rownames(conv.out),col.names=TRUE)
  }

  if(dic){
    cat("Start calculating deviance information criterion statistics...\n")
    dic.out<-dic.samples(model=jags.m,n.iter=n.iter,thin=n.thin)
    dev<-sum(dic.out$deviance)
    pen<-sum(dic.out$penalty)
    pen.dev<-dev+pen
    dic.stat<-rbind(dev,pen,pen.dev)
    rownames(dic.stat)<-c("D.bar","pD","DIC")
    colnames(dic.stat)<-""
    out$DIC<-dic.stat
  }

  if(mcmc.samples){
    out$mcmc.samples<-jags.out
  }

  if(!is.null(trace)){
    cat("Start saving trace plots...\n")
  }

  if(is.element("AR",trace)){
    for(i in 1:ntrt){
      png(paste("TracePlot_AR_",trtname[i],".png",sep=""),res=600,height=8.5,width=11,units="in")
      par(mfrow=c(n.chains,1))
      for(j in 1:n.chains){
        temp<-as.vector(jags.out[[j]][,paste("AR[",i,"]",sep="")])
        plot(temp,type="l",col="red",ylab="Absolute Risk",xlab="Iterations",main=paste("Chain",j))
      }
      dev.off()
    }
  }
  if(is.element("OR",trace)){
    for(i in 1:ntrt){
      for(k in 1:ntrt){
        if(i != k){
          png(paste("TracePlot_OR_",trtname[i],"_",trtname[k],".png",sep=""),res=600,height=8.5,width=11,units="in")
          par(mfrow=c(n.chains,1))
          for(j in 1:n.chains){
            temp<-as.vector(jags.out[[j]][,paste("OR[",i,",",k,"]",sep="")])
            plot(temp,type="l",col="red",ylab="Odds Ratio",xlab="Iterations",main=paste("Chain",j))
          }
          dev.off()
        }
      }
    }
  }
  if(is.element("LOR",trace)){
    for(i in 1:ntrt){
      for(k in 1:ntrt){
        if(i < k){
          png(paste("TracePlot_LOR_",trtname[i],"_",trtname[k],".png",sep=""),res=600,height=8.5,width=11,units="in")
          par(mfrow=c(n.chains,1))
          for(j in 1:n.chains){
            temp<-as.vector(jags.out[[j]][,paste("LOR[",i,",",k,"]",sep="")])
            plot(temp,type="l",col="red",ylab="Log Odds Ratio",xlab="Iterations",main=paste("Chain",j))
          }
          dev.off()
        }
      }
    }
  }
  if(is.element("RR",trace)){
    for(i in 1:ntrt){
      for(k in 1:ntrt){
        if(i != k){
          png(paste("TracePlot_RR_",trtname[i],"_",trtname[k],".png",sep=""),res=600,height=8.5,width=11,units="in")
          par(mfrow=c(n.chains,1))
          for(j in 1:n.chains){
            temp<-as.vector(jags.out[[j]][,paste("RR[",i,",",k,"]",sep="")])
            plot(temp,type="l",col="red",ylab="Risk Ratio",xlab="Iterations",main=paste("Chain",j))
          }
          dev.off()
        }
      }
    }
  }
  if(is.element("LRR",trace)){
    for(i in 1:ntrt){
      for(k in 1:ntrt){
        if(i < k){
          png(paste("TracePlot_LRR_",trtname[i],"_",trtname[k],".png",sep=""),res=600,height=8.5,width=11,units="in")
          par(mfrow=c(n.chains,1))
          for(j in 1:n.chains){
            temp<-as.vector(jags.out[[j]][,paste("LRR[",i,",",k,"]",sep="")])
            plot(temp,type="l",col="red",ylab="Log Risk Ratio",xlab="Iterations",main=paste("Chain",j))
          }
          dev.off()
        }
      }
    }
  }
  if(is.element("RD",trace)){
    for(i in 1:ntrt){
      for(k in 1:ntrt){
        if(i < k){
          png(paste("TracePlot_RD_",trtname[i],"_",trtname[k],".png",sep=""),res=600,height=8.5,width=11,units="in")
          par(mfrow=c(n.chains,1))
          for(j in 1:n.chains){
            temp<-as.vector(jags.out[[j]][,paste("RD[",i,",",k,"]",sep="")])
            plot(temp,type="l",col="red",ylab="Risk Difference",xlab="Iterations",main=paste("Chain",j))
          }
          dev.off()
        }
      }
    }
  }

  if(postdens){
    cat("Start saving posterior density plot for absolute risk...\n")
    mcmc<-NULL
    dens<-matrix(0,ntrt,3)
    colnames(dens)<-c("ymax","xmin","xmax")
    for(i in 1:ntrt){
      temp<-NULL
      for(j in 1:n.chains){
        temp<-c(temp,as.vector(jags.out[[j]][,paste("AR[",i,"]",sep="")]))
      }
      mcmc[[i]]<-temp
      tempdens<-density(temp)
      dens[i,]<-c(max(tempdens$y),quantile(temp,0.001),quantile(temp,0.999))
    }
    ymax<-max(dens[,"ymax"])
    xmin<-min(dens[,"xmin"])
    xmax<-max(dens[,"xmax"])
    cols<-rainbow(ntrt,s=1,v=0.6)
    pdf("AbsoluteRiskDensityPlot.pdf")
    par(mfrow=c(1,1),mar=c(5.5,5.5,2,2)+0.1)
    plot(density(mcmc[[1]]),xlim=c(xmin,xmax),ylim=c(0,ymax),xlab="Absolute Risk",ylab="Density",main="",col=cols[1],lty=1,lwd=2,cex.axis=2,cex.lab=2)
    for(i in 2:ntrt){
      lines(density(mcmc[[i]]),col=cols[i],lty=i,lwd=2)
    }
    legend("topright",legend=trtname,col=cols,lty=1:ntrt,lwd=2,cex=1.5)
    dev.off()
  }
  class(out)<-"nma.ab"
  return(out)
  options(warn=0)
}
