nma.ab.cont <-
function(s.id,t.id,mean,sd,total.n,data,trtname,param=c("mu","diff","rank.prob"),model="het_cor",prior.type,a=0.001,b=0.001,c=10,higher.better=FALSE,digits=4,n.adapt=5000,n.iter=100000,n.burnin=floor(n.iter/2),n.chains=3,n.thin=max(1,floor((n.iter-n.burnin)/100000)),conv.diag=FALSE,trace=NULL,dic=FALSE,postdens=FALSE,mcmc.samples=FALSE){
  ## check the input parameters
  options(warn=1)
  if(missing(s.id)) stop("need to specify study id.")
  if(missing(t.id)) stop("need to specify treatment.")
  if(missing(mean) | missing(sd)) stop("need to specify mean and sd of the continuous outcomes.")
  if(missing(total.n)) stop("need to specify total number.")
  if(!missing(data)){
    s.id<-eval(substitute(s.id),data,parent.frame())
    t.id<-eval(substitute(t.id),data,parent.frame())
    mean<-eval(substitute(mean),data,parent.frame())
    sd<-eval(substitute(sd),data,parent.frame())
    total.n<-eval(substitute(total.n),data,parent.frame())
  }
  if(length(s.id)!=length(t.id) | length(t.id)!=length(mean) | length(mean)!=length(sd) | length(sd)!=length(total.n) | length(total.n)!=length(s.id)){
    stop("s.id, t.id, mean, sd, and total.n have different lengths.")
  }
  if(!all(total.n>0)) stop("total number must be positive.")
  if(!all(total.n%%1==0)) warning("at least one event number or total number is not integer.")
  if(!is.element(model,c("hom_eqcor","het_eqcor","het_cor"))) stop("model should be specified as \"hom_eqcor\", \"het_eqcor\", or \"het_cor\".")

  if(any(is.na(s.id))|any(is.na(t.id))|any(is.na(mean))|any(is.na(sd))|any(is.na(total.n))){
    dat<-cbind(s.id,t.id,mean,sd,total.n)
    s.id<-s.id[complete.cases(dat)]
    t.id<-t.id[complete.cases(dat)]
    mean<-mean[complete.cases(dat)]
    sd<-sd[complete.cases(dat)]
    total.n<-total.n[complete.cases(dat)]
    cat("NA is not allowed in the input data set;\n")
    cat("the rows containing NA are removed.\n")
  }

  if(any(sd<=0)){
    s.id<-s.id[sd>0]
    t.id<-t.id[sd>0]
    mean<-mean[sd>0]
    total.n<-total.n[sd>0]
    sd<-sd[sd>0]
    cat("At least one sd is smaller than or equal to 0;\n")
    cat("the rows containing sd <= 0 are removed.\n")
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
    modelstring<-model.cont.hom.eqcor(prior.type,is.element("rank.prob",param))
  }
  if(model=="het_eqcor"){
    modelstring<-model.cont.het.eqcor(prior.type,is.element("rank.prob",param))
  }
  if(model=="het_cor"){
    modelstring<-model.cont.het.cor(prior.type,is.element("rank.prob",param))
  }

  ## jags data
  if(model == "hom_eqcor"| model == "het_eqcor"){
    if(prior.type == "unif"){
       if(is.element("rank.prob",param)) data.jags<-list(s=s.id,t=t.id,mean=mean,sd=sd,n=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),c=c,higher.better=higher.better)
       if(!is.element("rank.prob",param)) data.jags<-list(s=s.id,t=t.id,mean=mean,sd=sd,n=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),c=c)
    }
    if(prior.type == "invgamma"){
      if(is.element("rank.prob",param)) data.jags<-list(s=s.id,t=t.id,mean=mean,sd=sd,n=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),a=a,b=b,higher.better=higher.better)
      if(!is.element("rank.prob",param)) data.jags<-list(s=s.id,t=t.id,mean=mean,sd=sd,n=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),a=a,b=b)
    }
  }
  if(model=="het_cor"){
    if(prior.type == "invwishart"){
      I <- diag(ntrt)
      if(is.element("rank.prob",param)) data.jags<-list(s=s.id,t=t.id,mean=mean,sd=sd,n=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),I=I,higher.better=higher.better)
      if(!is.element("rank.prob",param)) data.jags<-list(s=s.id,t=t.id,mean=mean,sd=sd,n=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),I=I)
    }
    if(prior.type == "chol"){
      if(is.element("rank.prob",param)) data.jags<-list(s=s.id,t=t.id,mean=mean,sd=sd,n=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),c=c,higher.better=higher.better)
      if(!is.element("rank.prob",param)) data.jags<-list(s=s.id,t=t.id,mean=mean,sd=sd,n=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),c=c)
    }
  }

  ## jags initial value
  rng.seeds<-sample(1000000,n.chains)
  mu.init<-numeric(ntrt)
  for(i in 1:ntrt){
    mu.init[i]<-sum(mean[t.id==t.id[i]]*total.n[t.id==t.id[i]])/sum(total.n[t.id==t.id[i]])
  }
  init.jags<-list(NULL)
  if(model=="hom_eqcor"){
    if(prior.type=="unif"){
      for(ii in 1:n.chains){
        init.jags[[ii]]<-list(mu=mu.init,vi=matrix(0,nstudy,ntrt),sigma=c/2,rho=0.5,.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
      }
    }
    if(prior.type=="invgamma"){
      for(ii in 1:n.chains){
        init.jags[[ii]]<-list(mu=mu.init,vi=matrix(0,nstudy,ntrt),inv.sig.sq=a/b,rho=0.5,.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
      }
    }
  }
  if(model=="het_eqcor"){
    if(prior.type=="unif"){
      for(ii in 1:n.chains){
        init.jags[[ii]]<-list(mu=mu.init,vi=matrix(0,nstudy,ntrt),sigma=rep(c/2,ntrt),rho=0.5,.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
      }
    }
    if(prior.type=="invgamma"){
      for(ii in 1:n.chains){
        init.jags[[ii]]<-list(mu=mu.init,vi=matrix(0,nstudy,ntrt),inv.sig.sq=rep(a/b,ntrt),rho=0.5,.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
      }
    }
  }
  if(model=="het_cor"){
    if(prior.type=="invwishart"){
      for(ii in 1:n.chains){
        init.jags[[ii]]<-list(mu=mu.init,vi=matrix(0,nstudy,ntrt),T=(ntrt+1)*I,.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
      }
    }
    if(prior.type=="chol"){
      for(ii in 1:n.chains){
        init.jags[[ii]]<-list(mu=mu.init,vi=matrix(0,nstudy,ntrt),sigma=rep(c/2,ntrt),psi=matrix(3.1415926/2,ntrt-1,ntrt-1),.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
      }
    }
  }

  ## parameters to be monitored in jags
  if(!is.element("mu",param)) param<-c("mu",param)
  if(!is.null(trace)){
    if(!any(is.element(trace, param))) stop("at least one effect size in argument trace is not specified in argument param.")
  }
  monitor<-param[!is.element(param,c("diff"))]
  if(is.element("diff",param)){
    for(ii in 1:ntrt){
      for(jj in 1:ntrt){
        if(ii < jj) monitor<-c(monitor,paste("diff[",ii,",",jj,"]",sep=""))
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
  out$model<-"Normal likelihood with identity link."
  mu.id<-grep("mu",rownames(smry))
  mu.stat<-array(paste(format(round(smry[mu.id,"Mean"],digits=digits),nsmall=digits)," (",format(round(smry[mu.id,"SD"],digits=digits),nsmall=digits),")",sep=""),dim=c(ntrt,1))
  colnames(mu.stat)<-"Mean (SD)"
  rownames(mu.stat)<-trtname
  mu.quan<-array(paste(format(round(smry[mu.id,"50%"],digits=digits),nsmall=digits)," (",format(round(smry[mu.id,"2.5%"],digits=digits),nsmall=digits),", ",format(round(smry[mu.id,"97.5%"],digits=digits),nsmall=digits),")",sep=""),dim=c(ntrt,1))
  colnames(mu.quan)<-"Median (95% CI)"
  rownames(mu.quan)<-trtname
  out$TrtEffect<-list(Mean_SD=noquote(mu.stat),Median_CI=noquote(mu.quan))

  if(is.element("diff",param)){
    diff.stat<-diff.quan<-array(NA,dim=c(ntrt,ntrt))
    colnames(diff.stat)<-colnames(diff.quan)<-rownames(diff.stat)<-rownames(diff.quan)<-trtname
    for(i in 1:ntrt){
      diff.stat[i,i]<-diff.quan[i,i]<-"--"
      for(j in 1:ntrt){
        if(i < j){
          diff.ij<-paste("diff[",i,",",j,"]",sep="")
          diff.stat[i,j]<-paste(format(round(smry[diff.ij,"Mean"],digits=digits),nsmall=digits)," (",format(round(smry[diff.ij,"SD"],digits=digits),nsmall=digits),")",sep="")
          diff.stat[j,i]<-paste(format(round(-smry[diff.ij,"Mean"],digits=digits),nsmall=digits)," (",format(round(smry[diff.ij,"SD"],digits=digits),nsmall=digits),")",sep="")
          diff.quan[i,j]<-paste(format(round(smry[diff.ij,"50%"],digits=digits),nsmall=digits)," (",format(round(smry[diff.ij,"2.5%"],digits=digits),nsmall=digits),", ",format(round(smry[diff.ij,"97.5%"],digits=digits),nsmall=digits),")",sep="")
          diff.quan[j,i]<-paste(format(round(-smry[diff.ij,"50%"],digits=digits),nsmall=digits)," (",format(round(-smry[diff.ij,"97.5%"],digits=digits),nsmall=digits),", ",format(round(-smry[diff.ij,"2.5%"],digits=digits),nsmall=digits),")",sep="")
        }
      }
    }
    out$EffectDiff<-list(Mean_SD=noquote(diff.stat),Median_CI=noquote(diff.quan))
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

  if(is.element("mu",trace)){
    for(i in 1:ntrt){
      png(paste("TracePlot_mu_",trtname[i],".png",sep=""),res=600,height=8.5,width=11,units="in")
      par(mfrow=c(n.chains,1))
      for(j in 1:n.chains){
        temp<-as.vector(jags.out[[j]][,paste("mu[",i,"]",sep="")])
        plot(temp,type="l",col="red",ylab="Treatment Effect",xlab="Iterations",main=paste("Chain",j))
      }
      dev.off()
    }
  }
  if(is.element("diff",trace)){
    for(i in 1:ntrt){
      for(k in 1:ntrt){
        if(i < k){
          png(paste("TracePlot_diff_",trtname[i],"_",trtname[k],".png",sep=""),res=600,height=8.5,width=11,units="in")
          par(mfrow=c(n.chains,1))
          for(j in 1:n.chains){
            temp<-as.vector(jags.out[[j]][,paste("diff[",i,",",k,"]",sep="")])
            plot(temp,type="l",col="red",ylab="Effect Difference",xlab="Iterations",main=paste("Chain",j))
          }
          dev.off()
        }
      }
    }
  }

  if(postdens){
    cat("Start saving posterior density plot for treatment effect...\n")
    mcmc<-NULL
    dens<-matrix(0,ntrt,3)
    colnames(dens)<-c("ymax","xmin","xmax")
    for(i in 1:ntrt){
      temp<-NULL
      for(j in 1:n.chains){
        temp<-c(temp,as.vector(jags.out[[j]][,paste("mu[",i,"]",sep="")]))
      }
      mcmc[[i]]<-temp
      tempdens<-density(temp)
      dens[i,]<-c(max(tempdens$y),quantile(temp,0.001),quantile(temp,0.999))
    }
    ymax<-max(dens[,"ymax"])
    xmin<-min(dens[,"xmin"])
    xmax<-max(dens[,"xmax"])
    cols<-rainbow(ntrt,s=1,v=0.6)
    pdf("TreatmentEffectDensityPlot.pdf")
    par(mfrow=c(1,1),mar=c(5.5,5.5,2,2)+0.1)
    plot(density(mcmc[[1]]),xlim=c(xmin,xmax),ylim=c(0,ymax),xlab="Treatment Effect",ylab="Density",main="",col=cols[1],lty=1,lwd=2,cex.axis=2,cex.lab=2)
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