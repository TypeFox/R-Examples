GSA=function(x,y, genesets, genenames,
 method=c("maxmean","mean","absmean"), resp.type=c("Quantitative","Two class unpaired","Survival","Multiclass", "Two class paired", "tCorr", "taCorr"),
censoring.status=NULL,random.seed=NULL,  knn.neighbors=10,
  s0=NULL, s0.perc=NULL,minsize=15,maxsize=500, restand=TRUE,restand.basis=c("catalog","data"), nperms=200, xl.mode=c("regular","firsttime","next20","lasttime"), xl.time=NULL, xl.prevfit=NULL){
#
# computes feature set scores for a single set of data


this.call=match.call()
 method <- match.arg(method)
resp.type=match.arg(resp.type)
xl.mode=match.arg(xl.mode)
restand.basis=match.arg(restand.basis)

fdr.lo=NULL
fdr.hi=NULL
pvalues.lo=NULL
pvalues.hi=NULL

if(!is.null(random.seed)){
  set.seed(random.seed)
}
if(xl.mode=="regular" | xl.mode=="firsttime"){

if(sum(is.na(x))>0){
require(impute)
x=impute.knn(x,k=knn.neighbors)
}


#Error check: make sure that genenames and genesets have reasonable overlap


temp=match(unlist(genesets),genenames)
if(sum(!is.na(temp))/length(temp) < .05){
    stop("Fewer than 5% of genes in the genesets appear in the dataset. Make sure
that gene identifiers in dataset are Gene symbols")
}

junk=GSA.func(x,y,genesets=genesets, genenames=genenames,
             method=method, resp.type=resp.type,
 censoring.status=censoring.status,
              s0=s0,s0.perc=s0.perc, minsize=minsize,maxsize=maxsize,
           restand=restand, restand.basis= restand.basis)

r.obs=junk$score
stand.info=junk$stand.info
gene.scores=junk$gene.scores
s0=junk$s0
s0.perc=junk$s0.perc

r.star=NULL

k=length(genesets)
o=unlist(lapply(genesets,length))
gs.ind=(1:k)[o>= minsize &o<= maxsize]



catalog=NULL
ngenes=rep(NA,length(gs.ind))
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

# initialize for xl.iterations

r.star=matrix(0,nrow=length(genesets),ncol=nperms)
stand.info.star=matrix(0,nrow=8,ncol=nperms)
dimnames(stand.info.star)=list(c("mean.all","mean.abs","sd.all" , "sd.abs" , "mean.pos","sd.pos",  "mean.neg", "sd.neg"),NULL)

fdr.lo=NULL
fdr.hi=NULL
GSA.scores=NULL 
GSA.scores.perm=NULL
pvalues.lo=NULL
pvalues.hi=NULL
first.time=TRUE
}

if(xl.mode=="next20" |  xl.mode=="lasttime"){
 # get stuff from prevfit
GSA.scores=NULL
GSA.scores.perm=NULL
x=xl.prevfit$x
y=xl.prevfit$y
genesets=xl.prevfit$genesets
genenames=xl.prevfit$genenames
r.obs=xl.prevfit$r.obs
r.star=xl.prevfit$r.star
stand.info.star=xl.prevfit$stand.info.star
gs.mat=xl.prevfit$gs.mat
gs.ind=xl.prevfit$gs.ind
catalog=xl.prevfit$catalog
catalog.unique=xl.prevfit$catalog.unique
ngenes=xl.prevfit$ngenes
nperms=xl.prevfit$nperms
stand.info=xl.prevfit$stand.info
gene.scores=xl.prevfit$gene.scores
s0=xl.prevfit$s0
s0.perc=xl.prevfit$s0.perc
}



if(xl.mode=="regular"){
    first=1;last=nperms
  }
  if(xl.mode=="firsttime"){
    first=1;last=1
  }
  if(xl.mode=="next20"){
    first=xl.time; last= min(xl.time+9, nperms-1)
  }
  if(xl.mode=="lasttime"){
    first=nperms;last=nperms
  }




  for(i in first:last){
    if(i%%10==0){
    cat(c("perm=",i,paste("/",as.character(nperms)),sep=""),fill=T)
  }
     if(resp.type!="Two class paired"){oo=sample(1:length(y))}
     if(resp.type=="Two class paired"){oo=paired.perm(y)}

	 junk=GSA.func(x[,oo],y,genesets=genesets, genenames=genenames, 
              first.time=FALSE, return.gene.ind=FALSE, gs.mat=gs.mat, ngenes=ngenes,gs.ind=gs.ind, catalog=catalog,catalog.unique=catalog.unique,
             method=method, resp.type=resp.type, censoring.status=censoring.status,
               s0=s0, s0.perc=s0.perc,minsize=minsize,maxsize=maxsize,
           restand=restand)
      r.star[,i]=junk$score
      if (restand==TRUE) {  
          stand.info.star[, i] = unlist(junk$stand.info)
        } 
        else stand.info.star[, i] = NA   

  }

if(xl.mode=="regular" | xl.mode=="lasttime"){

 k=length(genesets)
pvalues.hi=rep(NA,k)
pvalues.lo=rep(NA,k)
for(i in gs.ind){
  pvalues.hi[i]=sum(r.star[i,]>r.obs[i])/nperms
  pvalues.lo[i]=sum(r.star[i,]<r.obs[i])/nperms
}

# estimate FDRs by plug-in method

cutp=c(0.001, 0.005,0.010,0.020,0.025,0.050,0.100,0.250,0.400,0.500)

res1=NULL;
for(i in 1:length(cutp)){
  o=!is.na(pvalues.lo)
  res1=rbind(res1,c(cutp[i], sum(pvalues.lo[o]<=cutp[i]), length(gs.ind)*cutp[i]))
}

f1=pmin(res1[,3]/res1[,2],1)
 
for(i in (length(f1)-1):1){
  f1[i]=min(f1[i],f1[i+1])
}

res2=NULL;
for(i in 1:length(cutp)){
 o=!is.na(pvalues.hi)
  res2=rbind(res2,c(cutp[i], sum(pvalues.hi[o]<=cutp[i]), length(gs.ind)*cutp[i]))
}

f2=pmin(res2[,3]/res2[,2],1)

for(i in (length(f2)-1):1){
  f2[i]=min(f2[i],f2[i+1])
}

fdr.lo=cbind(res1,f1)
fdr.hi=cbind(res2,f2)

fdr.lo=round(fdr.lo,3)
fdr.hi=round(fdr.hi,3)


dimnames(fdr.lo)=list(NULL, c("pv cutpoint","Number observed","Number expected","FDR"))
dimnames(fdr.hi)=list(NULL, c("pv cutpoint","Number observed","Number expected","FDR"))

GSA.scores=r.obs
GSA.scores.perm=r.star
# the last time through,  we delete  and other stuff that was needed just for
# the xl iters, 

x=NULL
y=NULL
gs.mat=NULL
genesets=NULL
genenames=NULL
catalog=NULL
catalog.unique=NULL
r.obs=NULL
r.star=NULL
}

out=list(
GSA.scores=GSA.scores,
GSA.scores.perm=GSA.scores.perm,
fdr.lo=fdr.lo,fdr.hi=fdr.hi, 
pvalues.lo=pvalues.lo,
pvalues.hi=pvalues.hi,
stand.info=stand.info,
stand.info.star=stand.info.star,
ngenes=ngenes,
nperms=nperms,
gene.scores=gene.scores,
s0=s0,
s0.perc=s0.perc,
resp.type=resp.type,
call=this.call,
x=x,
y=y,
genesets=genesets,
genenames=genenames,
r.obs=r.obs,
r.star=r.star,
gs.mat=gs.mat,
gs.ind=gs.ind,
catalog=catalog,
catalog.unique=catalog.unique)
class(out)="GSA"
out
}
