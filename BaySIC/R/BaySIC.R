#Libraries
library(fields)
library(rjags)
library(poibin)

#Character string DNA complementation# *DONE*
revcomp<-function(dna.seq){
  n.seq<-length(dna.seq)
  out.seq<-rep(NA,n.seq)
  for(i in 1:n.seq){
    out.seq[i]<-chartr("ACTG","TGAC",paste(rev(unlist(strsplit(dna.seq[i],""))),sep="",collapse=""))
    }
  return(out.seq)
  }

#Data Formatting#  
baysic.data<-function(dat,ref.dat,plot=FALSE,N=NULL,silent=TRUE){
  ref.check<-is.null(dim(ref.dat))
  if(ref.check==T){
    genes<-ref.dat[[1]][,1]
    seq.names<-colnames(ref.dat[[1]])[-1]
    }else{
    genes<-ref.dat[,1]
    seq.names<-colnames(ref.dat)[-c(1:2)]
    }
  n.g<-length(genes)
  if(silent==FALSE)
    dat<-dat[which(dat[,5]%in%c("silent","Silent","synonymous","Synonymous")==FALSE),]
  snv.dat<-matrix(0,ncol=32,nrow=n.g)
  indel.dat<-matrix(NA,ncol=1,nrow=n.g)
  dat.snv<-dat[dat[,5]!="INDEL",]
  dat.indel<-dat[dat[,5]=="INDEL",]
  
  #Case check
  case.check<-grep("[a|c|t|g]",as.character(dat.snv[,7]))
  if(length(case.check)>0)
    dat.snv[case.check,7]<-chartr("actgACTG","ACTGACTG",as.character(dat.snv[case.check,7]))
  seq.names<-chartr("actgACTG","ACTGACTG",seq.names)
  
  #Complementation check
  base.check<-grep("[A|G]",substr(as.character(dat.snv[,7]),2,2))
  if(length(base.check)>0)
    dat.snv[base.check,7]<-revcomp(as.character(dat.snv[base.check,7]))
  base.check.sn<-grep("[A|G]",substr(seq.names,2,2))
  if(length(base.check.sn)>0)
    seq.names[base.check.sn]<-revcomp(seq.names[base.check.sn])
      
  snv.fact<-factor(dat.snv[,7],levels=seq.names)
  
  #Adjust ref.dat colnames
  if(ref.check==T){
    for(i in 1:length(ref.dat)){
      colnames(ref.dat[[1]])[,-c(1:2)]<-seq.names
      }
    }else{
    colnames(ref.dat)[-c(1:2)]<-seq.names
    }
  if(length(unique(seq.names))<32)
    stop("Please check ref.dat, there are less than 32 unique trinucleotide motifs")
  for(i in 1:length(genes)){
    if(sum(as.character(dat.snv[,6])==as.character(genes[i]))>0)
      snv.dat[i,]<-table(snv.fact[which(as.character(dat.snv[,6])==as.character(genes[i]))])
    indel.dat[i]<-sum(as.character(dat.indel[,6])==as.character(genes[i]))
    }
  if(is.null(N)==T){
    if(ref.check==FALSE){
      N<-length(unique(dat[,4]))}else{
      N<-length(ref.dat)
      }
    }
  if(plot==TRUE){
    ids<-as.character(unique(dat[,4]))
    mut.type<-as.character(levels(dat[,5]))
    n.type<-length(mut.type)
    plot.dat<-matrix(0,nrow=n.type,ncol=length(ids))
    for(i in 1:length(ids)){
      plot.dat[,i]<-as.vector(table(dat[dat[,4]==ids[i],5]))
      }
    plot.dat<-plot.dat[,rev(order(apply(plot.dat,2,sum)))]
    layout(matrix(c(1,2,2,2),nrow=1))
    plot.dat2<-matrix(table(dat[,5])/nrow(dat),ncol=1)
    barplot(plot.dat2,border=NA,xlab="Mutations",ylab="Percentage of Mutations",col=rainbow(n.type))
    barplot(plot.dat,col=rainbow(n.type),legend.text=mut.type,ylab="Number of Mutations",space=0,xlab='Samples',border=NA)
    }
  colnames(snv.dat)<-seq.names
  out<-list("all.dat"=dat,"ref.dat"=ref.dat,"N"=N,"genes"=genes,"snv.dat"=snv.dat,"indel.dat"=indel.dat)
  return(out)
  }
  
##Plot mutation data##
BMR.plot<-function(dat.out){
  ref.dat<-dat.out$ref.dat
  snv.dat<-as.matrix(dat.out$snv.dat)
  N<-dat.out$N
  if(is.null(dim(ref.dat))){
    ref.dat<-as.matrix(Reduce("+",lapply(ref.dat,function(x) x[,-c(1:2)])))}else{
    ref.dat<-N*as.matrix(ref.dat[,-c(1:2)])
    }
  
  #Error Stops
  if(missing(snv.dat))
    stop("Please specify mutation data")
  if(missing(ref.dat))
    stop("Please specify reference data")
  if(nrow(snv.dat)!=nrow(ref.dat))
    stop("Dimensions of mutation data and reference data do not match")
  if(ncol(snv.dat)!=32|ncol(ref.dat)!=32)
    stop("Mutation and/or reference data is not in 32 sequence motif format")
  high.mut<-sum(snv.dat>ref.dat)
  if(high.mut>0)
    stop("Number of mutations greater than maximum possible per ref.dat")
  if(is.null(colnames(snv.dat))|is.null(colnames(ref.dat)))
    stop("Column names of data are not specified")

  mut.names<-colnames(snv.dat)
  ref.names<-colnames(ref.dat)
  if(sum(mut.names!=ref.names)>0)
    stop("Error:  Column names do not match")

  tri.revcomp<-revcomp(mut.names)
  CT.colname<-c()
  for(i in 1:32){
    if(substr(colnames(ref.dat)[i],2,2)%in%c("C","T")){
      CT.colname[i]<-colnames(ref.dat)[i]} else{
      CT.colname[i]<-tri.revcomp[i]}
    }
  out.mat<-matrix(NA,nrow=4,ncol=8)
  c.count<-apply(snv.dat,2,sum)/sum(snv.dat)
  ref.count<-apply(ref.dat,2,sum)/sum(as.numeric(ref.dat))
  count.ratio<-c.count/ref.count
  out.base<-c("A","C","T","G")
  in.base<-c("C","T")
  c.mat<-matrix(NA,4,4)
  t.mat<-matrix(NA,4,4)
  for(j in 1:4){
    for(k in 1:4){
      c.tri<-paste(out.base[j],"C",out.base[k],sep="")
      c.mat[j,k]<-count.ratio[which(CT.colname==c.tri)]
      t.tri<-paste(out.base[j],"T",out.base[k],sep="")
      t.mat[j,k]<-count.ratio[which(CT.colname==t.tri)]
      }
    }
  out.mat<-cbind(c.mat,t.mat)
  out.range<-range(log10(out.mat+.1))
  sig.colors<-rev(heat.colors(10))
  image.plot(x=1:ncol(out.mat),y=1:nrow(out.mat),col=sig.colors,z=log10(t(out.mat+.1)),xlab="",ylab="5' Base",axes=F,horizontal=T,legend.width=0.75,legend.shrink=0.50)
  title("")
  mtext("3' Base",side=3,line = par("mgp")[1])
  mtext(in.base,at=c(2.5,6.5),side=1,line = par("mgp")[2])
  mtext(out.base,at=1:4,side=2,line=par("mgp")[2],cex=0.75)
  mtext(rep(out.base,2),at=1:8,side=3,line=par("mgp")[2],cex=0.75)
  }


#Procedural JAGS model scripting *Done*
write.baysic<-function(mut.dat,covar=NULL,prior=NULL,fn.jags="baysic.jags"){
 SNV.cat<-colnames(mut.dat)[-ncol(mut.dat)]
 param.snv<-c("b.0",paste("b.",SNV.cat[-1],sep=""))
 param.string<-paste("b.0+",paste(param.snv[-1],"*X[i,",1:(length(SNV.cat)-1),"]",collapse="+",sep=""),sep="")
 if(is.null(covar)==FALSE){
  param.covar<-paste("g.",colnames(covar),sep="")
  param.snv<-c(param.snv,param.covar)
  param.string<-paste(param.string,"+",paste(param.covar,"*covar[gene[i],",1:ncol(covar),"]",collapse="+",sep=""),sep="")
  }
 model.string<-c("model{",
  "for(i in 1:n.mut){",
  paste("y[i]~dpois(exp(",param.string,")*N[i])}",sep="")
  ,"z~dpois(lambda*N.tot)")
 param<-c(param.snv,"lambda")
 if(is.null(prior)==FALSE){
  if(length(prior)!=length(param))
    stop("Length of prior argument does not match number of parameters")
	prior.string<-paste(param," ~ ",prior,sep="")
	}else{
	 prior.string<-c(paste(param[-length(param)]," ~ dnorm(0,0.001)",sep=""),"lambda ~ dgamma(0.001,0.001)")
	 }
 final.string<-c(model.string,"#Priors",prior.string,"}")
 fileconn<-file(fn.jags)
 writeLines(final.string,fileconn)
 close(fileconn)
 }


##SNV Category subroutine
fn.cat<-function(dat,snv.cat){
  if(is.list(snv.cat)==FALSE)
    stop("snv.cat improperly defined (not a list)")
  cat.check<-unlist(snv.cat)
  if(length(cat.check)!=32||sum(cat.check%in%1:32)!=32)
    stop("snv.cat improperly defined (assignments not mutually exclusive and/or not all assignments made") 
  n.cat<-length(snv.cat)
  if(is.null(names(snv.cat))==TRUE){
    cat.names<-paste("SNV_CAT_",1:n.cat,sep="")}else{
    cat.names<-names(snv.cat)
    }
  n.cat<-length(snv.cat)
  dat.cat<-matrix(NA,nrow=nrow(dat),ncol=n.cat)
  for(i in 1:n.cat){
    dat.cat[,i]<-apply(dat[,snv.cat[[i]]],1,sum)
    }
  colnames(dat.cat)<-cat.names
  return(dat.cat)
  }
    
## Model Fit ## *Done*
baysic.fit<-function(dat.out,snv.cat,covar=NULL,
  excl.list=NULL,burn.in=10000,n.samp=25000,fn.jags="baysic.jags",prior=NULL){
  genes<-dat.out$genes
  n.g<-length(genes)
  ref.dat<-dat.out$ref.dat
  if(is.null(dim(ref.dat))==TRUE){
    ref.dat<-as.matrix(Reduce("+",lapply(ref.dat,function(x) x[,-c(1:2)])))}else{
    N<-dat.out$N
    ref.dat<-N*as.matrix(ref.dat[,-c(1:2)])
    }
  snv.dat<-dat.out$snv.dat
  indel.dat<-dat.out$indel.dat
  
  #Redefine data through SNV categories
  snv.dat.cat<-fn.cat(snv.dat,snv.cat)
  ref.dat.cat<-fn.cat(ref.dat,snv.cat)
  n.cat<-ncol(snv.dat.cat)
  
  #Remove excluded genes (if any)
  if(is.null(excl.list)==FALSE){
    if(is.numeric(excl.list)==TRUE){
      snv.dat.cat<-snv.dat.cat[-excl.list,]
      ref.dat.cat<-ref.dat.cat[-excl.list,]
      }
    else{
      which.excl<-which(genes%in%excl.list)
      snv.dat.cat<-snv.dat.cat[-which.excl,]
      ref.dat.cat<-ref.dat.cat[-which.excl,]
      }
    }
  
  mut.dat<-cbind(snv.dat.cat,indel.dat)
  
  #Write BaySIC model
  write.baysic(mut.dat,covar,prior,fn.jags)
  
  #Define BaySIC Data
  if(is.null(covar)==FALSE){
    X<-rbind(matrix(0,nrow=n.g,ncol=n.cat-1),kronecker(diag(rep(1,n.cat-1)),matrix(1,nrow=n.g,ncol=1)))
    dat.jags<-list("n.mut"=nrow(X),"X"=X,"gene"=rep(1:nrow(snv.dat.cat),n.cat),"covar"=covar,"y"=as.vector(snv.dat.cat),"z"=sum(indel.dat),"N"=as.vector(ref.dat.cat),"N.tot"=sum(as.numeric(ref.dat.cat)))
    }else{
    X<-rbind(rep(0,n.cat-1),diag(rep(1,n.cat-1)))
    dat.jags<-list("n.mut"=nrow(X),"X"=X,"y"=apply(snv.dat.cat,2,sum),"z"=sum(indel.dat),"N"=apply(ref.dat.cat,2,function(x) sum(as.numeric(x))),"N.tot"=sum(as.numeric(ref.dat.cat)))
    }
    
  #Fit BaySIC model
  baysic.model<-jags.model(fn.jags,
    data = dat.jags,
    n.chains = 2,
    n.adapt = 1000)
  update(baysic.model,burn.in)
  jags.param<-c("b.0",paste("b.",colnames(snv.dat.cat)[-1],sep=""),"lambda")
  if(is.null(covar)==FALSE)
    jags.param<-c(jags.param,paste("g.",colnames(covar)))
  baysic.post<-coda.samples(baysic.model,jags.param,n.samp)
  out<-list("fit.post"=baysic.post,"covar"=covar,"snv.cat"=snv.cat,"excl.list"=excl.list)
  return(out)
  }

#Fuzzy FDR# *Done*  
fuzzy.FDR.approx<-function(pprev,p,alpha,N){
  g<-length(p)
  out.mat<-matrix(NA,nrow=g,ncol=N)
  for(i in 1:g){if(pprev[i]>p[i]){pprev[i]<-p[i]}}
  for(i in 1:N){
    out.mat[,i]<-p.adjust(runif(g,pprev,p),"fdr")
    }
  out.prob<-apply(out.mat,1,function(x) sum(x<=alpha)/N)
  return(out.prob)
  }

#Baysic Testing#
baysic.test<-function(dat.out,fit.out,fdr.level=0.15,fuzzy.cnt=10000,r=NULL,subtype=NULL,PB.approx=FALSE){ 
  all.dat<-dat.out$all.dat
  N<-dat.out$N
  genes<-dat.out$genes
  n.g<-length(genes)
  ref.dat<-dat.out$ref.dat
  ref.check<-is.null(dim(ref.dat))
  covar<-fit.out$covar                 
  snv.cat<-fit.out$snv.cat
  n.cat<-length(snv.cat)
  
  baysic.post<-as.matrix(fit.out$fit.post)
  n.param.samp<-nrow(baysic.post)
  if(is.null(r)==FALSE)
    baysic.post<-baysic.post[seq.int(1,n.param.samp,length.out=r),]
  which.lambda<-which(colnames(baysic.post)=="lambda")
  post.snv<-baysic.post[,-which.lambda]
  post.lambda<-baysic.post[,which.lambda]
  if(PB.approx==TRUE){
    PB.method<-"RNA"}else{
    PB.method<-"DFT-CF"
    }
  if(is.null(subtype)==FALSE&ref.check==TRUE){
    sub.ids<-subtype[,1]
    ref.ids<-names(ref.dat)
    if(sum(sub.ids%in%ref.ids)!=nrow(subtype))
      stop("Ids in ref.dat list do not match ids in subtype data")
    }
  n.mut<-rep(NA,n.g)
  for(i in 1:length(genes)){
    n.mut[i]<-length(unique(all.dat[as.character(all.dat[,6])==as.character(genes[i]),4]))
    }
  if(ref.check==TRUE){
    ref.dat.cat<-lapply(ref.dat,function(x) fn.cat(x[,-c(1:2)],snv.cat=snv.cat))
    ref.dat.chr<-ref.dat[[1]][,2]}else{
    ref.dat.cat<-fn.cat(ref.dat[,-c(1:2)],snv.cat)
    ref.dat.chr<-ref.dat[,2]
    }
  
  #Subtype Data
  if(is.null(subtype)==FALSE){
    all.dat<-dat.out$all.dat
    sub.names<-unique(subtype[,2])
    n.sub<-length(sub.names)
    sub.nmuts<-matrix(NA,nrow=n.g,ncol=n.sub)
    sub.res<-matrix(NA,nrow=n.g,ncol=3*n.sub)
    sub.N<-rep(NA,n.sub)
    if(ref.check==TRUE){
      sub.ref.dat<-list()}else{
      sub.ref.dat<-ref.dat.cat
      }
    for(j in 1:n.sub){
      sub.N[j]<-sum(subtype[,2]==sub.names[j])
      which.j<-as.character(subtype[which(subtype[,2]==sub.names[j]),1])
      dat.j<-all.dat[which(all.dat[,4]%in%which.j),]
      for(i in 1:n.g){
        sub.nmuts[i,j]<-length(unique(dat.j[as.character(dat.j[,6])==as.character(genes[i]),4]))
        }
      if(ref.check==TRUE)
        sub.ref.dat[[j]]<-ref.dat.cat[which(ref.ids%in%which.j)]
      }
    }
    
  p.out<-matrix(NA,nrow=n.g,ncol=2)
  if(is.null(covar)==FALSE){
    for(i in 1:n.g){
      cov.post.i<-post.snv[,-c(1:(n.cat))]%*%t(covar[i,,drop=F])
      prob.snv.i<-exp(cbind(post.snv[,1],sweep(post.snv[,2:n.cat,drop=F],1,post.snv[,1,drop=F]+cov.post.i,"+")))
      T.x.i<-n.mut[i]
      rate.i<-cbind(prob.snv.i,post.lambda)
      if(ref.check==TRUE){
        ref.i<-sapply(ref.dat.cat,function(x) c(x[i,],sum(x[i,])))
        p.i<-rate.i%*%ref.i
        p.mut<-1-ppois(0,p.i)
        p.out[i,1]<-ifelse(T.x.i>0,mean(1-apply(p.mut,1,function(x) ppoibin(T.x.i-1,x,method=PB.method))),1)
        p.out[i,2]<-mean(1-apply(p.mut,1,function(x) ppoibin(T.x.i,x,method=PB.method)))
        if(is.null(subtype)==FALSE){
          for(j in 1:n.sub){
            T.x.ij<-sub.nmuts[i,j]
            N.j<-sub.N[j]
            ref.ij<-sapply(sub.ref.dat[[j]],function(x) c(x[i,],sum(x[i,])))
            p.ij<-rate.i%*%ref.ij
            p.mut.j<-1-ppois(0,p.ij)
            sub.res[i,3*(j-1)+1]<-ifelse(T.x.ij>0,mean(1-apply(p.mut.j,1,function(x) ppoibin(T.x.ij-1,x,method=PB.method))),1)
            sub.res[i,3*(j-1)+2]<-mean(1-apply(p.mut.j,1,function(x) ppoibin(T.x.ij,x,method=PB.method)))
            }
          }
        }else{
        ref.i<-matrix(c(ref.dat.cat[i,],sum(ref.dat.cat[i,])),ncol=1)
        p.i<-rate.i%*%ref.i
        p.mut<-1-ppois(0,p.i)
        p.out[i,1]<-ifelse(T.x.i>0,mean(1-pbinom(T.x.i-1,N,p.mut)),1)
        p.out[i,2]<-mean(1-pbinom(T.x.i,N,p.mut))
        if(is.null(subtype)==FALSE){
          for(j in 1:n.sub){
            T.x.ij<-sub.nmuts[i,j]
            N.j<-sub.N[j]
            ref.ij<-matrix(c(sub.ref.dat[i,],sum(sub.ref.dat[i,])),ncol=1)
            p.ij<-rate.i%*%ref.ij
            p.mut.j<-1-ppois(0,p.ij)
            sub.res[i,3*(j-1)+1]<-ifelse(T.x.ij>0,mean(1-pbinom(T.x.ij-1,N.j,p.mut.j)),1)
            sub.res[i,3*(j-1)+2]<-mean(1-pbinom(T.x.ij,N.j,p.mut.j))
            }
          }
        }
      }
    }else{
    prob.snv<-exp(cbind(post.snv[,1],sweep(post.snv[,-1,drop=F],1,post.snv[,1,drop=F],"+")))
    rate.all<-cbind(prob.snv,post.lambda)
    for(i in 1:n.g){
      T.x.i<-n.mut[i]
      if(ref.check==TRUE){
        ref.i<-sapply(ref.dat.cat,function(x) c(x[i,],sum(x[i,])))
        p.i<-rate.all%*%ref.i
        p.mut<-1-ppois(0,p.i)
        p.out[i,1]<-ifelse(T.x.i>0,mean(1-apply(p.mut,1,function(x) ppoibin(T.x.i-1,x,method=PB.method))),1)
        p.out[i,2]<-mean(1-apply(p.mut,1,function(x) ppoibin(T.x.i,x,method=PB.method)))
        if(is.null(subtype)==FALSE){
          for(j in 1:n.sub){
            T.x.ij<-sub.nmuts[i,j]
            N.j<-sub.N[j]
            ref.ij<-sapply(sub.ref.dat[[j]],function(x) c(x[i,],sum(x[i,])))
            p.ij<-rate.all%*%ref.ij
            p.mut.j<-1-ppois(0,p.ij)
            sub.res[i,3*(j-1)+1]<-ifelse(T.x.ij>0,mean(1-apply(p.mut.j,1,function(x) ppoibin(T.x.ij-1,x,method=PB.method))),1)
            sub.res[i,3*(j-1)+2]<-mean(1-apply(p.mut.j,1,function(x) ppoibin(T.x.ij,x,method=PB.method)))
            }
          }
        }else{
        ref.i<-matrix(c(ref.dat.cat[i,],sum(ref.dat.cat[i,])),ncol=1)
        p.i<-rate.all%*%ref.i
        p.mut<-1-ppois(0,p.i)
        p.out[i,1]<-ifelse(T.x.i>0,mean(1-pbinom(T.x.i-1,N,p.mut)),1)
        p.out[i,2]<-mean(1-pbinom(T.x.i,N,p.mut))
        if(is.null(subtype)==FALSE){
          for(j in 1:n.sub){
            T.x.ij<-sub.nmuts[i,j]
            N.j<-sub.N[j]
            ref.ij<-matrix(c(sub.ref.dat[i,],sum(sub.ref.dat[i,])),ncol=1)
            p.ij<-rate.all%*%ref.ij
            p.mut.j<-1-ppois(0,p.ij)
            sub.res[i,3*(j-1)+1]<-ifelse(T.x.ij>0,mean(1-pbinom(T.x.ij-1,N.j,p.mut.j)),1)
            sub.res[i,3*(j-1)+2]<-mean(1-pbinom(T.x.ij,N.j,p.mut.j))
            }
          }
        }
      } 
    }
  #Obtain fuzzy FDR estimates
  fuzzy.res<-fuzzy.FDR.approx(p.out[,2],p.out[,1],fdr.level,fuzzy.cnt)
  test.res<-data.frame(as.character(genes),as.character(ref.dat.chr),n.mut,p.out,fuzzy.res)
  colnames(test.res)<-c("gene","chr","N_mut","pppval","pppval+","fuzzy.prob")
  if(is.null(subtype)==FALSE){
    colnames(sub.res)<-rep("X",n.sub*3)
    for(j in 1:n.sub){
      sub.res[,3*j]<-fuzzy.FDR.approx(sub.res[,3*(j-1)+2],sub.res[,3*(j-1)+1],fdr.level,fuzzy.cnt)
      colnames(sub.res)[(3*(j-1)+1):(3*j)]<-paste(c("pppval","pppval+","fuzzy.prob"),sub.names[j],sep="_")
      }
    test.res<-cbind(test.res,sub.res)
    }
  return(list(test.res=test.res,fdr.level=fdr.level,fuzzy.cnt=fuzzy.cnt,subtype=subtype))
  }

  
