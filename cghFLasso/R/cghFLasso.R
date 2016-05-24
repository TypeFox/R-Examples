#################################################################
############################## basic functions
#################################################################

mystand <- function(x)
{
        # standardizes columns so that  mean 0, var= 1/n (as in lars)
        nm <- dim(x)
        n <- nm[1]
        m <- nm[2]
        im <- inactive <- seq(m)
        one <- rep(1, n)
        vn <- dimnames(x)[[2]]
        ### Center x and y, and scale x, and save the means and sds
        meanx <- drop(one %*% x)/n
        x <- scale(x, meanx, FALSE)
        # centers x
        normx <- sqrt(drop(one %*% (x^2)))
        names(normx) <- NULL
        x <- scale(x, FALSE, normx)
        return(x=x, meanx=meanx, normx=normx)
}

###############
countdf <- function(beta)
{
        nval <- 0
        if(beta[1] != 0) {
                nval <- 1
        }
        for(i in 2:length(beta)) {
                if((beta[i] != beta[i - 1]) & beta[i] != 0) {
                        nval <- nval + 1
                }
        }
        return(nval)
}

##################
sumabs <- function(beta){
   sum(abs(beta))
 }

###################
sumabsd <- function(beta){

        abs(beta[1]) + sum(abs(diff(beta)))
}



###################
my.mad<-function(x.v)
{
x.m<-median(x.v)
result<-median(abs(x.v-x.m))
return(result)
}


##################################################################
###############################  Fused Lasso by Nike
##################################################################

L2L1Vit = function(obsSeq, lambda2 = 1, lambda1 = 0, maxSegs = 1000,
                   retPrimBds = F)
{    
  if(!is.double(obsSeq)){
    obsSeq = as.double(obsSeq)
  }
  obsWts = rep(1.0, length(obsSeq))                            
  nSegs = rep( as.integer(0) , length(obsSeq))
  midSegs = matrix(nrow=2, ncol=length(obsSeq), data=0)
  
  if(!is.double(lambda2)) lambda2 = as.double(lambda2)
  if(!is.double(lambda1)) lambda1 = as.double(lambda1)
  
  if(retPrimBds){
    primBds = rep(0.0, 2)
  }else{
    primBds = NULL
  }
  
  retBeta = rep(0, length(obsSeq))
    
  .Call("L2L1Vit", 
          obsSeq, obsWts, lambda2, lambda1, 
          retBeta, maxSegs, nSegs, 
          midSegs, primBds)
          
  return(list(beta=retBeta, vitMsgs=NULL, backPtrs=midSegs, nVitMsgSegs=nSegs,
              lambda2 = lambda2, lambda1 = lambda1, primBds = primBds) )
}

####################################
L2L1VitPath = function(obsSeq, lambda2 = c(1.0,3.0), lambda1 = 0, maxSegs = 1000,
                       segmentedFit = T)
{ 
  if(!segmentedFit) warning("L2L1VitPath : segmentedFit = F is not efficient.")
  
  if(!is.double(obsSeq)){
    obsSeq = as.double(obsSeq)
  }
  if(!is.double(lambda2)){
    lambda2 = as.double(lambda2)
  }          
  if(!is.double(lambda1)) lambda1 = as.double(lambda1)
  
  primBds = matrix(nrow=2, ncol=length(lambda2), data=0)
  
  if(segmentedFit){
    betaPath = vector(mode="list", length=length(lambda2))
    betaSegs = vector(mode="list", length=length(lambda2))
  }else{
    betaPath = matrix(nrow=length(obsSeq), ncol=length(lambda2), data=0)
    betaSegs = NULL
  } 
  
  .Call("L2L1VitPath", obsSeq, lambda2, lambda1, 
        betaPath, maxSegs, betaSegs, primBds)
  return(list(betaPath=betaPath, betaSegs=betaSegs, lambda2 = lambda2, 
              lambda1 = lambda1,
              primBds = primBds)) 
}

######################################
L2L1VitPathPrimBds = function(retList, lambda1, subFit = NULL)
{
  nFit = length( retList$betaPath )
  
  if(is.null(subFit) ){
    subFit = 0:(nFit-1)
  }else{
    subFit = as.integer(subFit - 1)
    if(max(subFit) >= nFit || min(subFit) < 0) stop("index out of range.");
  }
  
  if(!is.double(lambda1)) lambda1 = as.double(lambda1)
  
  lambda1 = lambda1 - retList$lambda1
  
  if(min(lambda1) < 0) stop("lambda1 must be non-negative")
  
  retArr = array(0.0, c(2, length(lambda1), length(subFit))  )
  
  .Call("L2L1GetPrimBds", 
        retList$betaPath, retList$betaSegs, 
        lambda1, subFit, retArr)
        
  return( retArr )
}

###############################################
L2L1ExpandFit = function(retList, lambda1 = 0.0, listInds = 1)
{   
  if( length(lambda1) > 1){
    warning("only first value of lambda1 will be used")
  }
    
  if(is.matrix(retList$betaSegs) ){
     retList = L2L1PathExtract(retList, listInds)
     return( L2L1ExpandFit( retList  = retList,
                            lambda1  = lambda1,
                            listInds = 1:length(listInds) ) ) 
  }   
  nFit = length( retList$betaPath )
  
  seqLen = retList$betaSegs[[1]]
  seqLen = seqLen[2, ncol(seqLen)]
  
  listInds = as.integer(listInds - 1)
  if(max(listInds) >= nFit || min(listInds) < 0) stop("index out of range.");

  if( length(listInds) == 1){
    retBeta = rep(0.0, seqLen)
  }else{
    retBeta = matrix( nrow=seqLen, ncol=length(listInds), data=0)
  }   
  
  lambda1 = lambda1 - retList$lambda1
  
  .Call("L2L1ExpandFit", 
        retList$betaPath, retList$betaSegs, lambda1, 
        listInds, retBeta)
  return( retBeta )
}


######################################
MatchPrimBd = function(o, lam1seq, lam2seq, s1, s2){

  # Get the path
  flpth = L2L1VitPath(o, lambda2 = lam2seq, lambda1 = 0.0, segmentedFit = T)

  # Check the bounds for different values of lambda1
  primBds = L2L1VitPathPrimBds(flpth, lambda1 = lam1seq, subFit = NULL)
  
  # Find the nearest match :
  err1 = ( abs( primBds[1,,] - s1 ) / s1 )
  err2 = ( abs( primBds[2,,] - s2 ) / s2 )

  i1 = which.min(err1 + err2) - 1
  lam2i = 1 + as.integer( i1 / length(lam1seq) )
  lam1i = 1 + ( i1 %% length(lam1seq) )

  # Get the final fit for that value of lambda1/lambda2

  bv1 = L2L1ExpandFit(flpth, lambda1 = lam1seq[lam1i], listInds = lam2i)

  return( list(beta    = bv1, 
               s1      = primBds[1,lam1i,lam2i], 
               s2      = primBds[2,lam1i,lam2i],
               lam1    = lam1seq[lam1i], 
               lam2    = lam2seq[lam2i],
               lam1i   = lam1i,
               lam2i   = lam2i,
               err1    = err1,
               err2    = err2,
               flpth   = flpth,
               primBds = primBds ))
}


####################################################################
############################## application on CGH
####################################################################

###/////////////////////////////////////////////////////////////////
########## For calculating region FDR: test mean=0 for each segments 

segment.score<-function(y.v, beta.v)
{ 
  ###### y.v: one chromosome of the original CGH array
  ###### beta.v: fused.lasso result of this chromosome
  
  gap<-diff(beta.v)
  gap.point<-c(0, (1:length(gap))[gap!=0], length(y.v))
  k<-length(gap.point)
  result<-numeric(k-1)
  size<-numeric(k-1)
  sum.var<-0
  for(i in 2:k)
    {
     begin<-gap.point[i-1]+1  
     end<-gap.point[i]
     cur<-y.v[begin:end]
     size[i-1]<-end-begin+1
     result[i-1]<-sum(beta.v[begin:end])/size[i-1]^0.5
     sum.var<-sum.var+sum((cur-mean(cur))^2)
    }
   emp.sd<-(sum.var/length(y.v))^0.5
   result<-result/emp.sd
   pvalue<-1-pnorm(abs(result))
   return(list(pvalue=pvalue, score=result, emp.sd=emp.sd, size=size))
}

############ Q-value for a given p-value cutoff
fdr.one<-function(cut, seg.obj)
{
  obs.count<-sum(seg.obj$size[seg.obj$pvalue<=cut])
  null.count<-sum(seg.obj$size)*cut
  result<-null.count/obs.count
}

########### Q-value vector for each spot of the array
fdr.seg<-function(seg.obj)
{
   temp<-seg.obj$pvalue
   qvalue<-apply(as.matrix(temp), 1, fdr.one, seg.obj=seg.obj)
   result<-rep(qvalue, seg.obj$size)
   return(result)
}

####//////////////////////////////////////////////////////////////////////////////
########## function to call alteration using fused.constant.lasso on one chromosome

CGH.fused.lasso.one.fast<-function(y.v)
{
my.y<-y.v[!is.na(y.v)]

### for s2
y.lo<-lowess(my.y, f=min(1/20, 10/length(my.y)), delta=0)$y
temp.d<-diff(y.lo)
sd.diff<-my.mad(temp.d)
s2=2*sd.diff+sum(abs(temp.d)[abs(temp.d)>4*sd.diff])

### for s1
y.lo.s1<-lowess(my.y, f=max(1/50, 50/length(my.y)), delta=0)$y
s1=sum(abs(y.lo.s1))


################## get sequence
lam1.max=0.01
lam1.list=seq(0,lam1.max,by=0.001)
np1<-length(lam1.list)

lam2.max=max((length(my.y)/10)^0.5, 5)
dlm=0.05
lam2.min=min(0.05*(length(my.y)/10)^0.5, 1)
lam2.list=seq(lam2.min,lam2.max, by=dlm)

result1=MatchPrimBd(my.y, lam1.list, lam2.list, s1, s2)
beta<-result1$beta

result<-y.v
result[!is.na(y.v)]<-beta
return(result)
}


################# function to call alteration on one CGH array at a given FDR level

CGH.FusedLasso.One<-function(cgh, chromosome, FL.norm=NULL, FDR=NULL)
{
############ cgh: vector of CGH array log_2 ratio, ordered according to genome position
############ chromosome: vector of chromosome number of each gene in "cgh"
############ FL.norm:    vector of fused lasso result on normal reference array


Good<-!is.na(cgh)

cgh.y<-cgh[Good]
chrom<-chromosome[Good]
FL.beta<-cgh.y

ch.set<-sort(unique(chrom))


####### perform fused lasso on each chromosome
for(ch in ch.set)
 {
   print(ch)
   FL.beta[chrom==ch]<-CGH.fused.lasso.one.fast(cgh.y[chrom==ch])
 }

###### if no FDR specified 
if(is.null(FDR))
  {
   result<-rep(NA, length(cgh))
   result[Good]<-FL.beta
   return(result)
  }

###### otherwise, only return alteration survives the specified FDR level

if(!is.null(FL.norm))
 {
   #### control FDR based on reference arrays
   theta<-seq(0.01, 2, by=0.01)
   fdr.v<-theta
   for(i in 1:length(theta))
    {
     obs.count<-sum(abs(FL.beta)>=theta[i], na.rm=T)/sum(!is.na(FL.beta))
     exp.count<-sum(abs(FL.norm)>=theta[i], na.rm=T)/sum(!is.na(FL.norm))
     fdr.v[i]<-exp.count/(obs.count+0.0001)
    }

   FL.y<-FL.beta
   FL.y[abs(FL.beta)<min(theta[fdr.v<=FDR])]<-0
 } else {
   #### control FDR based on segments
    Qvalue<-rep(99, length(cgh.y))
    for(ch in ch.set)
     {
        cur.beta<-FL.beta[chrom==ch]
        cur.y<-cgh.y[chrom==ch]
        cur.seg<-segment.score(cur.y, cur.beta)
        cur.q<-fdr.seg(cur.seg)
        Qvalue[chrom==ch]<-cur.q
     }  
    FL.y<-FL.beta
    FL.y[Qvalue>=FDR]<-0
 }

result<-rep(NA, length(cgh))
result[Good]<-FL.y
return(result)
}

#####################

plugin.NA<-function(result, ch, size=5)
{
  WHOLEdata<-as.matrix(result)
  chset<-unique(ch)
  WHOLEdata[is.na(WHOLEdata)]<-999
  WHOLEdata.sm<-NULL

  p<-ncol(WHOLEdata)
  storage.mode(p) <- "integer"
  storage.mode(size)<-"integer"
 
#if(!is.loaded("avesmooth")){
#  dyn.load("CGH.FusedLasso.so")
#}  

  for(j in sort(chset)){
   ntemp<-sum(ch==j)
   if(ntemp>1)
     {
      temp<-WHOLEdata[ch==j, ]
      result<-temp
      storage.mode(ntemp)<-"integer"
      storage.mode(temp) <- "double" 
      storage.mode(result) <- "double"

      junk <- .Fortran("avesmooth",
			ntemp,
			p,
                        size, temp, result=result, 
                        PACKAGE="cghFLasso")

      y<-junk$result
      }
    if(ntemp==1)
      y<-WHOLEdata[ch==j,]

     WHOLEdata.sm<-rbind(WHOLEdata.sm, y)
  }
WHOLEdata.sm[WHOLEdata.sm==999]=0
result.sm<-WHOLEdata
result.sm[WHOLEdata==999]<-WHOLEdata.sm[WHOLEdata==999]
return(result.sm)
}


##################################################################
############################ function for users
#################################################################

cghFLasso.ref<-function(CGH.Array, chromosome=NULL, filter=NULL)
{
   CGH.Array<-as.matrix(CGH.Array)

   ### take a subset of the array
   if(is.null(filter))
      filter<-rep(0, nrow(CGH.Array))

   if(is.null(chromosome))
    {
        chromosome=rep(1, sum(filter==0))
    } else {
        chrom<-chromosome[filter==0]
    }
    
    
    array<-CGH.Array[filter==0,]
    chrom<-chromosome[filter==0]

    result<-apply(array, 2, CGH.FusedLasso.One, chromosome=chrom)
    result<-as.matrix(result) 
    return(result)
}

##############################################################
#############################  key function


cghFLasso<-function(CGH.Array, chromosome=NULL, nucleotide.position=NULL, FL.norm=NULL, FDR=NULL, filter=NULL, missing.PlugIn=TRUE, smooth.size=5)
{
   CGH.Array<-as.matrix(CGH.Array)

   ### take a subset of the array
   if(is.null(filter))
      filter<-rep(0, nrow(CGH.Array))

   if(is.null(chromosome))
    {
        chromosome=rep(1, sum(filter==0))
    } else {
        chromosome<-chromosome[filter==0]
    }

   if(is.null(nucleotide.position))
    {
        nucposi=1:sum(filter==0)
    } else {
        nucposi<-nucleotide.position[filter==0]
    }
    
    Array<-as.matrix(CGH.Array[filter==0,])

    result<-apply(Array, 2, CGH.FusedLasso.One, chromosome=chromosome, FL.norm=FL.norm, FDR=FDR)
    result<-as.matrix(result)

    ##### put in NA through average smooth
    if(missing.PlugIn)
    {
    result<-plugin.NA(result, chromosome, size=smooth.size)
    }

    colnames(result)<-colnames(as.matrix(CGH.Array))
 
    ans<-list(Esti.CopyN=result, CGH.Array=Array, 
              chromosome=chromosome, nucleotide.position=nucposi, FDR=FDR)
    class(ans)<-"cghFLasso"
    return(ans)
}

######################################################################################
################ function to summary
###############################################

consensusFDR<-function(fdr, ConsensusCount, s.n)
{
N<-length(ConsensusCount)
m<-s.n
p<-fdr
ccset<-sort(unique(ConsensusCount))

tempresult<-1:(max(ccset)+1)
for(i in ccset)
{
curcc<-i
nu<-1-pbinom(curcc-1, size=m, prob=p)
den<-sum(ConsensusCount>=curcc)/N
tempresult[curcc+1]<-nu/den
}
tempresult[ConsensusCount+1]
}

########
summary.cghFLasso<-function(object, index, ...)
{
 CGH.FL.obj=object
 s.n=length(index)
 if(s.n==1)
    stop("Summary should be applied to a group of arrays!")

 raw.Array<-CGH.FL.obj$CGH.Array[,index] 
 Esti.copy<-CGH.FL.obj$Esti.CopyN[,index] 

 ConsensusCount<-apply(Esti.copy!=0, 1, sum, na.rm=T)
 CC.FDR=consensusFDR(CGH.FL.obj$FDR, ConsensusCount, s.n)

 Amp.CC=apply(Esti.copy>0, 1, sum, na.rm=T)
 Del.CC=apply(Esti.copy<0, 1, sum, na.rm=T)  

 chromosome=CGH.FL.obj$chromosome
 nucposi=CGH.FL.obj$nucleotide.position
 
 ch.set=unique(chromosome)
 chrom.summary<-matrix(NA, nrow=length(ch.set), ncol=4)
 for(i in 1:length(ch.set))
  {
    cur.ch<-ch.set[i]
    cur.pick<-chromosome==cur.ch
    chrom.summary[i,1]=cur.ch
    chrom.summary[i,2]=max(Amp.CC[cur.pick])/s.n
    chrom.summary[i,3]=max(Del.CC[cur.pick])/s.n
    chrom.summary[i,4]=mean(CC.FDR[cur.pick]<0.05)
  }
 colnames(chrom.summary)<-c("Chrom", "Max.Amp.Sample%", "Max.Del.Sample%", "Alter.Chromosome%")

 sample.summary<-apply(Esti.copy!=0, 2, mean, na.rm=T)

 ans<-list(ConsensusCount=ConsensusCount, CC.FDR=CC.FDR, 
           Amp.ConCount=Amp.CC, Del.ConCount=Del.CC, 
           chrom.summary=chrom.summary, sample.summary=sample.summary,
           Esti.copy=Esti.copy, chromosome=chromosome, nucposi=nucposi)
 class(ans)<-"summary.cghFLasso"
 ans
}

##########
print.summary.cghFLasso<-function(x, ...)
{
  x.obj=x
  chrom.summary=x.obj$chrom.summary
  pick<-apply(chrom.summary[,2:4]>0, 1, sum, na.rm=T)
  if(sum(pick>0)>1)
   {
    print(round(chrom.summary[pick>0,],3))
   } 
  if(sum(pick>0)==1) 
   {
      temp<-matrix(chrom.summary[pick>0,], nrow=1)
      colnames(temp)<-colnames(chrom.summary) 
      print(round(temp,3))
   }
  if(sum(pick>0)==0)
    {
      print("No Alteration!")      
    }
}

##########

output.cghFLasso<-function(summary.obj, file, gene.info=NULL)
{
  summary.obj<-unclass(summary.obj)
  temp<-summary.obj$Esti.copy 
  if(is.null(colnames(temp)))
    colnames(temp)<-paste("Sample", 1:ncol(temp))
  if(is.null(gene.info))
    gene.info=paste("Gene", 1:nrow(temp))
  output<-cbind(gene.info, summary.obj$chromosome,summary.obj$nucposi,summary.obj$ConsensusCount, summary.obj$CC.FDR, temp)
  output.2<-rbind( c(rep("", ncol(as.matrix(gene.info))), rep("", 4), round(summary.obj$sample.summary, 3)), output)
  output.2[1,1]<-"Alter.%"
  colnames(output.2)<-c(paste("Gene.Info", colnames(gene.info), sep="."), "Chromosome Number", "Nuc Position", "Consensus Count", "Consensus FDR", colnames(temp))
  write.table(output.2, file=file, row.names=FALSE, col.names=TRUE, sep="\t")
}



######################################################################################
################# function to plot 
###############################################


plot.CGH.FL.Single<-function(raw.Array, Esti.copy, chromosome, nucposi, centro)
{
  sample=raw.Array 
  clustindex=Esti.copy 
  chr=chromosome
  centr=centro
  graylevel=0.9

    n=length(chr)
    Ch=max(chr)
    ran=range(nucposi)
    plot(nucposi,rep(1,length(nucposi)),type="n",axes=F,ylim=c(0,max(chr)+1),
    xlim=c(ran[1],ran[2]),ylab="",xlab="")

upcol<-gray(0.5)
downcol<-gray(0.5)

con<-log(10)
sample<-sample/con *log(2)

chset<-unique(chr)

###########
## First, draw background scales
for(j in c(sort(chset), Ch+1)){
   jp=Ch-j+1
   nuc<-nucposi[chr==j | chr==j-1]
   xlim=range(0,max(nuc))
   for(copy in 2:10)
   segments(xlim[1], jp+log(copy)/con, xlim[2], jp+log(copy)/con, col=gray(graylevel))
}
##########
## Second, draw data
for(j in sort(chset)){
   jp=Ch-j+1
   nuc=nucposi[chr==j]
   y=sample[chr==j]
   indexy<-clustindex[chr==j]

   y[is.na(y)]<-0
   y[y==999]<-0
   yposi<-y
   yposi[y<0]<-0
   ynega<-y
   ynega[y>0]<-0

   pick<-(1:length(indexy))[indexy!=0 & (indexy*y>0)]
   npick<-(1:length(indexy))[indexy==0 | (indexy*y<=0)]

  ### plot unpicked one
   if(length(npick)>0)
    {
    segments(nuc[npick],jp,nuc[npick],jp+yposi[npick],col=upcol)
    segments(nuc[npick],jp,nuc[npick],jp+ynega[npick],col=downcol)
    }
  ### plot picked one
   if(length(pick)>0)
     {
    segments(nuc[pick],jp,nuc[pick],jp+yposi[pick],col=2)
    segments(nuc[pick],jp,nuc[pick],jp+ynega[pick],col=3)
     }
   xlim=range(0,max(nuc))
   segments(xlim[1],jp,xlim[2],jp)
   if(j<23)
      text(-5000000,jp,labels=j,cex=.7)
   if(j==23)
      text(-5000000, jp, labels="X", cex=0.7)
   if(j==24)
      text(-5000000, jp, labels="Y", cex=0.7)
   segments(centr[j],jp+0.25,centr[j],jp-0.25,col=6)
   }
}

##########################################################################
##########################################################################
plot.CGH.FL.Consensus<-function(raw.Array, Esti.copy, chromosome, nucposi, centro, cut=0)
{
sampleM=Esti.copy
chr=chromosome 
centr=centro
graylevel=0.9 

  ####
    n<-length(chr)
    ran<-range(nucposi)
 
   M<-ncol(sampleM)    
   color<-rainbow(6*M)
   height<-M   
   chset<-unique(chr)    

####### for legend location, if there are more than 20 chromosome(most probably from human genome), put legend in the middle
####### otherwise, put the legend at the bottom of the picture    
if(length(chset)>20)
{
plot(nucposi,rep(1,length(nucposi)),type="n",axes=F,ylim=c(0,max(chr)*2+1),
    xlim=c(ran[1],ran[2]),ylab="",xlab="")
}
if(length(chset)<20)
{
plot(nucposi,rep(1,length(nucposi)),type="n",axes=F,ylim=c(-7,max(chr)*2+1),
    xlim=c(ran[1],ran[2]),ylab="",xlab="")
}

Ch=max(chr)  
############
## first plot background scale lines
for(j in sort(chset)){
   jp=Ch*2-j*2+1
   nuc<-nucposi[chr==j]
   xlim=range(0,max(nuc))
   for(copy in seq(0.2, 0.8, 0.2))
   {
   segments(xlim[1], jp+copy, xlim[2], jp+copy, col=gray(graylevel))
   segments(xlim[1], jp-copy, xlim[2], jp-copy, col=gray(graylevel))
   }
   segments(xlim[1], jp+1, xlim[2], jp+1, col=gray(graylevel-0.1))
   segments(xlim[1], jp-1, xlim[2], jp-1, col=gray(graylevel-0.1))
}
############
## second plot data
for(j in sort(chset))
   {
   jp=Ch*2-j*2+1
   nuc=nucposi[chr==j]
   xlim=range(0,max(nuc))
   y=sampleM[chr==j, ]
   y[is.na(y)]=0  
   posicount<-apply(y>0, 1, sum)+1
   negacount<-apply(y<0, 1, sum)+1
   for(i in 1:length(nuc)){
      if(posicount[i]-1+negacount[i]-1>=cut)
        {
         segments(nuc[i],jp,nuc[i],jp+(posicount[i]-1)/height,col=color[M-posicount[i]])
      ######if(negacount[i]-1 >=cut)
         segments(nuc[i],jp,nuc[i],jp-(negacount[i]-1)/height,col=color[M+negacount[i]])
         }
    }
   if(j<23)
     text(-5000000,jp,labels=j,cex=.7)
   if(j==23)
     text(-5000000, jp, labels="X", cex=0.7)
   if(j==24)
     text(-5000000, jp, labels="Y", cex=0.7)
   segments(centr[j],jp+.5,centr[j],jp-.5,col=6)
   segments(xlim[1],jp,xlim[2],jp)                  
   }

###############
## Third. Plot legend

ends<-range(nucposi)
start<-ends[1]+(ends[2]-ends[1])*1/2
finish<-ends[1]+(ends[2]-ends[1])*6/7

if(length(chset)>20)
  level<-5*2
if(length(chset)<20)
  level<-(-4)

segments(start, level-2, start,level+2)
segments(start, level-2, finish, level-2)
segments(finish,level-2, finish, level+2)
segments(finish,level+2, start, level+2)

interval<-(finish-start)/40
M<-100
color<-rainbow(6*M)
height<-M   
xpoints<-seq(start+interval, finish-interval, length=2*M)

for(h in seq(0.2, 1, 0.2))
  {
   segments(xpoints[1], level+h, xpoints[2*M], level+h, col=gray(0.9))
   segments(xpoints[1], level-h, xpoints[2*M], level-h, col=gray(0.9))
  }
for(i in 1:M)
 {
  segments(xpoints[i+M], level, xpoints[i+M],level+i/height, col=color[M-i])
  segments(xpoints[M+1-i], level, xpoints[M+1-i],level-i/height, col=color[M+i])
 } 
segments(xpoints[1], level, xpoints[2*M], level, col=1)
gap<-M%/%2
yh<-level-3
text(xpoints[1], yh, 100)
text(xpoints[gap], yh, 50)
text((xpoints[M]+xpoints[M+1])/2, yh, 0)
text(xpoints[M+gap], yh, 50)
text(xpoints[2*M], yh, 100)
text(start-100, level+3, "Percent of Samples")

}

##########################################################################
##########################################################################
plot.CGH.FL.All<-function(raw.Array, Esti.copy, chromosome, nucposi)
{
  arrayraw=raw.Array
  arrayregion=Esti.copy
  chr=chromosome
  
  
  arrayname<-colnames(arrayraw)

  n<-length(chr)
  m<-ncol(arrayraw)
  plot(1:n,rep(1,n),type="n",axes=F,ylim=c(0,2*m+1), xlim=c(1,n),ylab="",xlab="")


upcol<-"#FFD600"
downcol<-"#EBFF00"
chcol<-"#00FFFF"

    for(j in 1:m)
     {
    jp<-2*m+1-2*j
   segments(0, jp+0.5,n, jp+0.5, col=upcol)
   segments(0, jp+1,n, jp+1, col=upcol)
   segments(0, jp-0.5,n, jp-0.5, col=downcol)
    y<-arrayraw[,j]
    y[y==999]<-0
    y[is.na(y)]<-0
    segments((1:n)[y>0],jp,(1:n)[y>0],jp+y[y>0],col=upcol)
    segments((1:n)[y<0],jp,(1:n)[y<0],jp+y[y<0],col=downcol)
     
    interval<-arrayregion[,j]
    pick<-(1:n)[interval>0]
    if(sum(pick)>0)
       segments((1:n)[pick], jp, (1:n)[pick], jp+y[pick], col=2)    
    pick<-(1:n)[interval<0]
    if(sum(pick)>0)
       segments((1:n)[pick], jp, (1:n)[pick], jp+y[pick], col=3)    

   text(-n/50,jp,labels=j,cex=.7)
   text(n/50, jp+1, labels=arrayname[j], cex=0.5)
   text(n/50, jp+0.5,labels=paste("Gain/Loss=", round(mean(arrayregion[,j]!=0),3)*100, "%", sep=""), cex=0.8)
     }

   chset<-unique(chr)
   count<-chset
   for(j in sort(chset))
      count[j]<-sum(chr==j)
   count<-cumsum(count)
   segments(count, -0, count, 2*m+1, col=chcol)
   count<-c(0, count)
   for(j in sort(chset))
      text((count[j]+count[j+1])/2, 2*m+1, j, cex=0.5, col=chcol)
}

#################
#################
plot.cghFLasso<-function(x, index, type="All", centro=NULL, ...)
{
CGH.FL.obj=x

if(type=="Single" | type=="Lines")
{
   cur.index<-index[1]
} else {
    cur.index=index
}

raw.Array<-CGH.FL.obj$CGH.Array[,cur.index] 
Esti.Copy<-CGH.FL.obj$Esti.CopyN[,cur.index] 

if(type=="Lines")
{
   plot(raw.Array, xlab="genome order",pch=20, col=gray(0.5), ylab="log_2 ratio")
   abline(h=0, col=gray(0.3))
   points(Esti.Copy, col=2, type="l", lwd=2)    
   return(type)
} 

raw.Array[is.na(raw.Array)]<-0
Esti.Copy[is.na(Esti.Copy)]<-0


chromosome=CGH.FL.obj$chromosome
nucposi=CGH.FL.obj$nucleotide.position


if(type=="All")
 {
  plot.CGH.FL.All(as.matrix(raw.Array), as.matrix(Esti.Copy), chromosome, nucposi)
  return(type)
 }


if(is.null(centro))
{
  centro<-c(123000000, 93500000, 91500000, 51000000, 47700000, 60500000, 58900000, 45000000, 49000000, 40000000,53000000, 35400000, 15000000, 15600000, 17000000, 39000000, 24000000, 16000000, 28500000, 27800000, 12200000, 11900000, 58500000, 10000000)
}

if(type=="Single")
  {
   plot.CGH.FL.Single(raw.Array, Esti.Copy, chromosome, nucposi, centro)
   return(type)
  }

if(type=="Consensus")
  {
   plot.CGH.FL.Consensus(raw.Array, Esti.Copy, chromosome, nucposi, centro)
   return(type)
  }

return(FALSE)
}



