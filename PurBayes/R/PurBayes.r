#PurBayes
library(rjags)

write.PB<-function(fn.jags,prior=NULL,het=FALSE,germ=FALSE){
  if(germ==TRUE){
    germ.string<-c("Z~dbin(q,M)",
      "q~dunif(0,1)")}else{
    germ.string<-""}
  if(is.null(prior)==TRUE)
    prior<-"dunif(0,1)"
  if(het==FALSE){ 
    model.string<-c("model{",germ.string,
      "for(i in 1:N.snv){",
      ifelse(germ==TRUE,"n.alt[i] ~ dbin(q*pur,n.tot[i])}","n.alt[i] ~ dbin(0.5*pur,n.tot[i])}"),
      paste("pur ~ ",prior,"}",sep=""))
      }else{
    model.string<-c("model{",
      germ.string,
      "for(j in 1:n.pop){",
      "lambda[j]~dunif(0,1)}",
      "lambda.srt<-sort(lambda)",
      "for(i in 1:N.snv){",
      "n.cat[i] ~ dcat(kappa)",
      ifelse(germ==TRUE,"n.alt[i] ~ dbin(q*lambda.srt[n.cat[i]],n.tot[i])}",
        "n.alt[i] ~ dbin(0.5*lambda.srt[n.cat[i]],n.tot[i])}"),
      "kappa ~ ddirich(alpha[])",
      "pur<-max(lambda)",
      "}"
      )
      }
	fileconn<-file(fn.jags)
  writeLines(model.string,fileconn)
  close(fileconn)
  }
  
plot.PurBayes<-function(x,...){
  plot(x$N,x$Y,cex=0.75,pch=16,xlab="Total Reads",ylab="Mutant Allele Reads",...)
  n.pop<-x$n.pop
  PB.post<-as.matrix(x$PB.post)
  N.max<-max(x$N)
  if(n.pop==1){
    which.pur<-which(colnames(PB.post)=="pur")
    val.j<-quantile(PB.post[,which.pur],c(0.025,0.5,0.975))
    lines(c(0,N.max),c(0,N.max*0.5*val.j[2]))
    lines(c(0,N.max),c(0,N.max*0.5*val.j[1]),lty=2)
    lines(c(0,N.max),c(0,N.max*0.5*val.j[3]),lty=2)
    }else{
    for(i in 1:n.pop){
      which.lambda.i<-which(colnames(PB.post)==paste("lambda.srt[",i,"]",sep=""))
      val.j<-quantile(PB.post[,which.lambda.i],c(0.025,0.5,0.975))
      lines(c(0,N.max),c(0,N.max*0.5*val.j[2]))
      lines(c(0,N.max),c(0,N.max*0.5*val.j[1]),lty=2)
      lines(c(0,N.max),c(0,N.max*0.5*val.j[3]),lty=2)
      }
    }
  }

dic.run<-function(pb.model){
  #First pass
  dic.popt<-dic.samples(pb.model,1000,type="popt")
  if(is.na(sum(dic.popt$penalty))==TRUE){
    print("Warning: Observation(s) with NaN popt values. Using 2*pD approximation")
    dic.pD<-dic.samples(pb.model,1000,type="pD")
    nan.list<-which(is.na(dic.popt$penalty))
    dic.popt$penalty[nan.list]<-dic.pD$penalty[nan.list]*2
    }
  return(dic.popt)
  }
  
PurBayes<-function(N,Y,M=NULL,Z=NULL,pop.max=5,prior=NULL,burn.in=50000,n.post=10000,fn.jags="PB.jags",plot=FALSE){
  n.pop<-1
  germ.dat<-ifelse(is.null(M)==F&&is.null(Z)==F,TRUE,FALSE)
  write.PB(fn.jags,prior,het=FALSE,germ=germ.dat)
  if(germ.dat==TRUE){
    pb.dat.old<-list("N.snv"=length(N),"n.tot"=N,"n.alt"=Y,"M"=sum(M),"Z"=sum(Z))}else{
    pb.dat.old<-list("N.snv"=length(N),"n.tot"=N,"n.alt"=Y)}
  
  pb.m.old<-jags.model(file=fn.jags,pb.dat.old,n.chains=2,n.adapt=1000)
  update(pb.m.old,burn.in)
  
  dic.old<-dic.run(pb.m.old)
  
  pb.list<-list()
  dic.list<-list()
  pb.list[[1]]<-pb.m.old
  dic.list[[1]]<-dic.old
  
  pop.max.break<-0
  repeat{
    if(n.pop==pop.max){
      print("Warning: PurBayes has reached the defined pop.max, consider increasing this value")
      pop.max.break<-1
      break
      }
    n.pop<-n.pop+1
    write.PB(fn.jags,prior,het=TRUE,germ=germ.dat)
    if(germ.dat==TRUE){
      pb.dat<-list("N.snv"=length(N),"n.tot"=N,"n.alt"=Y,"M"=sum(M),"Z"=sum(Z),"n.pop"=n.pop,"alpha"=rep(1,n.pop))}else{
      pb.dat<-list("N.snv"=length(N),"n.tot"=N,"n.alt"=Y,"n.pop"=n.pop,"alpha"=rep(1,n.pop))}
    pb.m<-jags.model(file=fn.jags,pb.dat,n.chains=2,n.adapt=1000)
    update(pb.m,burn.in)
    dic.new<-dic.run(pb.m)
    pb.list[[n.pop]]<-pb.m
    dic.list[[n.pop]]<-dic.new
    dic.check<-diffdic(dic.old,dic.new)
    dic.diff<-sum(dic.check)
  
    if(dic.diff<0)
      break
    pb.m.old<-pb.m
    dic.old<-dic.new
    }
  
  which.ref<-ifelse(pop.max.break==0,n.pop-1,n.pop)
  dic.ref<-dic.list[[which.ref]]
  
  #Identify optimal model  
  model.check<-matrix(NA,nrow=n.pop,ncol=4)
  model.check[which.ref,]<-c(which.ref,sum(dic.ref$deviance)+sum(dic.ref$penalty),0,0)
  for(i in c(1:n.pop)[-which.ref]){
    dic.i<-dic.list[[i]]  
    dic.check.i<-diffdic(dic.i,dic.ref)
    dic.diff.i<-sum(dic.check.i)
    dic.sd.i<-sqrt(length(dic.check.i))*sd(dic.check.i)
    model.check[i,]<-c(i,sum(dic.i$deviance)+sum(dic.i$penalty),dic.diff.i,dic.sd.i)
    }
  n.pop.fin<-min(which(model.check[1:(which.ref-1),3]<=model.check[1:(which.ref-1),4]),which.ref)
     
  colnames(model.check)<-c("N.pop","PED","Delta","SD(Delta)")
  print(paste("PurBayes detected ",n.pop.fin," population(s) of variants",sep=""))
  pb.m.final<-pb.list[[n.pop.fin]]
  if(n.pop.fin==1){
    if(germ.dat==FALSE){
      pb.post<-coda.samples(pb.m.final,"pur",n.post)}else{
      pb.post<-coda.samples(pb.m.final,c("pur","q"),n.post)}
    }else{
    if(germ.dat==FALSE){
      pb.post<-coda.samples(pb.m.final,c("pur","kappa","lambda.srt"),n.post)}else{
      pb.post<-coda.samples(pb.m.final,c("pur","kappa","lambda.srt","q"),n.post)}
    }
  names(pb.list)<-paste("n_",1:length(pb.list),sep="")
  out<-list("N"=N,"Y"=Y,"M"=M,"Z"=Z,"n.pop"=n.pop.fin,"dev.mat"=as.data.frame(model.check),"PB.post"=pb.post,"which.ref"=which.ref,"jags.fits"=pb.list)
  class(out)<-"PurBayes"
  if(plot==TRUE)
    plot.PurBayes(out)
  return(out)
  }
  
summary.PurBayes<-function(object,...){
  post.out<-t(apply(as.matrix(object$PB.post),2,function(x) quantile(x,c(0.50,0.025,0.975))))
  colnames(post.out)<-c("Median","2.5% Quant.","97.5% Quant.")
  out.dev<-object$dev.mat
  out.dev$Ref_Model<-rep("",nrow(out.dev))
  out.dev$Ref_Model[object$which.ref]<-"*"
  pur.dist<-post.out[which(row.names(post.out)=="pur"),]
  n.pop<-object$n.pop
  out<-list("purity"=pur.dist,"post.dist"=post.out,"n.pop"=n.pop,"dev.out"=out.dev)
  class(out)<-"summary.PurBayes"
  print(out)
  invisible(out)
  }
  
print.summary.PurBayes<-function(x,...){
  cat("Purity Estimate\n")
  print(x$purity)
  cat("\n")
  cat("Number of Populations\n")
  print(x$n.pop)
  cat("\n")
  cat("Model Posterior Distributions\n")
  print(x$post.dist)
  cat("\n")
  cat("Penalized Deviance Results\n")
  print(x$dev.out,row.names=F)
  }