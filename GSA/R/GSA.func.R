GSA.func=function(x,y, genesets, genenames, geneset.names=NULL,
 method=c("maxmean","mean","absmean"),
resp.type=c("Quantitative","Two class unpaired","Survival","Multiclass", "Two class paired", "tCorr", "taCorr"),
censoring.status=NULL,
first.time=TRUE,  return.gene.ind=TRUE, ngenes=NULL, gs.mat=NULL, gs.ind=NULL,
catalog=NULL, catalog.unique=NULL, s0=NULL, s0.perc=NULL, minsize=15,maxsize=500, restand=TRUE, restand.basis=c("catalog","data")){

#
# computes gene set scores for a single set of data

this.call=match.call()
 method <- match.arg(method)
 resp.type <- match.arg(resp.type)
restand.basis <-  match.arg(restand.basis)

BIG=10e9
me=rowMeans(x)

if(resp.type=="tCorr" | resp.type=="taCorr"){s0=0}


if(first.time){
k=length(genesets)
o=unlist(lapply(genesets,length))
gs.ind=(1:k)[o>= minsize &o<= maxsize]

ngenes=rep(NA,length(gs.ind))

# note that gs.mat is initialized with the value nrow(x)+1; this is
# a trick: see note below
catalog=NULL
gs.mat=matrix(nrow(x)+1,nrow=length(gs.ind),ncol=maxsize)
ii=0
for(i in gs.ind){
  ii=ii+1
  gene.set=match(genesets[[i]],genenames)
  gene.set=gene.set[!is.na(gene.set)]
 catalog=c(catalog,gene.set)
 if(length(gene.set)>0){
 gs.mat[ii,1:length(gene.set)]=gene.set
 ngenes[ii]=length(gene.set)
 }
}

catalog.unique=unique(catalog)
 gs.mat=gs.mat[,1:max(ngenes, na.rm=TRUE)]
}


# estimate s0 if necessary

  if(is.null(s0)){
    
     if(!is.null(s0.perc)){
      if((s0.perc != -1 & s0.perc < 0) | s0.perc > 100){
       stop("Illegal value for s0.perc: must be between 0 and 100, or equal
to (-1) (meaning that s0 should be set to zero)")
       }

     if(s0.perc== -1){s0=0}
    }
     
initflag=FALSE
if(is.null(s0.perc)){initflag=TRUE}
else if(s0.perc>=0){initflag=TRUE}

 if(initflag){
       
  if(resp.type=="Quantitative"){
    init.fit=quantitative.func(x[catalog.unique,],y,s0=0)
      }
  if(resp.type=="Two class unpaired"){
    init.fit=ttest.func(x[catalog.unique,],y,s0=0)
      }
   if(resp.type=="Survival"){
   init.fit=cox.func(x[catalog.unique,],y,censoring.status,s0=0)
 }
    if(resp.type=="Multiclass"){
     init.fit=multiclass.func(x[catalog.unique,],y,s0=0)
  } 
 if(resp.type=="Two class paired"){
    init.fit=paired.ttest.func(x[catalog.unique,],y,s0=0)
      }

}
     

     if(is.null(s0.perc)){
      s0=est.s0(init.fit$tt,init.fit$sd)$s0.hat
       s0.perc=100*sum(init.fit$sd<s0)/length(init.fit$sd)
     }
     else if (s0.perc!=-1){ s0 <- quantile(init.fit$sd,s0.perc/100)}
   }


###  end of estimation of s0

if(resp.type=="Two class unpaired"){
  jun=ttest.func(x[catalog.unique,],y,s0=s0)
  if(restand & restand.basis=="data"){
     jun.stand=ttest.func(x,y,s0=s0)
  }
}



if(resp.type=="Survival"){
 jun=cox.func(x[catalog.unique,],y,censoring.status,s0=s0)
 if(restand & restand.basis=="data"){
     jun.stand=cox.func(x,y,s0=s0)
  }
}

if(resp.type=="Multiclass"){
 jun=multiclass.func(x[catalog.unique,],y,s0=s0)
  if(restand & restand.basis=="data"){
     jun.stand=multiclass.func(x,y,s0=s0)
  }
}

if(resp.type=="Quantitative"){
  jun=quantitative.func(x[catalog.unique,],y,s0=s0)
  if(restand & restand.basis=="data"){
     jun.stand=quantitative.func(x,y,s0=s0)
  }
}
if(resp.type=="Two class paired"){
  jun=paired.ttest.func(x[catalog.unique,],y,s0=s0)
  if(restand & restand.basis=="data"){
     jun.stand=paired.ttest.func(x,y,s0=s0)
  }
}
if(resp.type=="tCorr"){
  jun=tCorr.func(x[catalog.unique,],y,s0=s0)
  if(restand & restand.basis=="data"){
   jun.stand=tCorr.func(x,y,s0=s0)
  }
}
if(resp.type=="taCorr"){
  jun=taCorr.func(x[catalog.unique,],y,s0=s0)
  if(restand & restand.basis=="data"){
   jun.stand=taCorr.func(x,y,s0=s0)
  }
}

if(restand.basis=="catalog"){
	tt=rep(NA,nrow(x))
 	 s=tt
	tt[catalog.unique]=jun$tt
	s[catalog.unique]=jun$sd
}

if(restand.basis=="data"){
  	tt = jun.stand$tt  
    	s = jun.stand$sd
}

gene.scores=tt


mean.all=mean.abs=sd.all=sd.abs=mean.pos=sd.pos=mean.neg=sd.neg=NULL

if(restand){
 if(restand.basis=="catalog"){
	mean.all=mean(tt[catalog])
	sd.all=sqrt(var(tt[catalog]))
	
	mean.abs=mean(abs(tt[catalog]))
	sd.abs=sqrt(var(abs(tt[catalog])))
	
	mean.pos=mean(tt[catalog]*(tt[catalog]>0))
	sd.pos=sqrt(var(tt[catalog]*(tt[catalog]>0)))
	
	mean.neg=-mean(tt[catalog]*(tt[catalog]<0))
	sd.neg=sqrt(var(tt[catalog]*(tt[catalog]<0)))
}
if(restand.basis=="data"){
     mean.all = mean(tt)
        sd.all = sqrt(var(tt))
        mean.abs = mean(abs(tt))
        sd.abs = sqrt(var(abs(tt)))
        mean.pos = mean(tt * (tt > 0))
        sd.pos = sqrt(var(tt * (tt > 0)))
        mean.neg = -mean(tt * (tt < 0))
        sd.neg = sqrt(var(tt * (tt < 0)))
    }

}



stand.info=list(mean.all=mean.all,
mean.abs=mean.abs,
sd.all=sd.all,
sd.abs=sd.abs,
mean.pos=mean.pos,
sd.pos=sd.pos,
mean.neg=mean.neg,
sd.neg=sd.neg
)

# Note: this is a trick: I artificially tack a zero on the end of tt,
# so that  tt2[gs.mat] has the same dim as gs.mat, with zeroes filling
# in the end of each row

tt2=c(tt,0)
ttt=matrix(tt2[gs.mat],nrow=nrow(gs.mat))



if(method=="maxmean"){
#rpos=rowSums(pmax(ttt,0))/ngenes
#rneg=-1*rowSums(pmin(ttt,0))/ngenes

s2=abs(ttt)
rpos=rowSums((ttt+s2)/2)/ngenes
rneg=rowSums((s2-ttt)/2)/ngenes

rpos[is.na(rpos)]=0
rneg[is.na(rneg)]=0
if(restand ){rpos=(rpos- mean.pos)/(sd.pos)}
if(restand & resp.type!= "Multiclass" & resp.type!= "taCorr"){rneg=(rneg- mean.neg)/(sd.neg)}
rr=pmax(rpos,rneg)
rr[rneg>rpos]=-1*rr[rneg>rpos]
}

gene.ind=NULL






if(method=="mean"){
 rr=rowSums(ttt)/ngenes
 if(restand){rr=(rr-mean.all)/(sd.all/sqrt(ngenes))}
}

if(method=="absmean"){
rr=rowSums(abs(ttt))/ngenes
 if(restand){rr=(rr-mean.all)/(sd.all/sqrt(ngenes))}
}

if(return.gene.ind){
 gene.ind=vector("list",length(gs.ind))
ii=0
for(i in gs.ind){
ii=ii+1
  gene.set=match(genesets[[i]],genenames)
  gene.set=gene.set[!is.na(gene.set)]

  if(rr[ii]>0){ gene.ind[[i]]=(1:length(gene.set))[tt[gene.set]>0]}
  if(rr[ii]<0){ gene.ind[[i]]=(1:length(gene.set))[tt[gene.set]<0]}
}}



rrr=rep(NA,length(genesets))
rrr[gs.ind]=rr
rrt=-qnorm(1-pt(rrr, df=length(y)-2))

out=list(scores=rrr, norm.scores=rrt,  mean=me, sd=s, gene.ind=gene.ind, geneset.names=geneset.names, gene.scores=gene.scores,s0=s0,s0.perc=s0.perc,
stand.info=stand.info,
method=method,call=this.call)
class(out)="GSA.func"
out
}


   

