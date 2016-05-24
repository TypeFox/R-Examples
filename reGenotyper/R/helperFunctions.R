.deltaTmatAtOneMK<-function(trait,gt1,t0.mki,fileName="test",indSample=NULL) {
  #this function is used by the function .deltaTMatAllmk_list(see below)
  #using matrixeQTL do not need consider all gt is 1 or 0 case, it will still gives t value ~0 and p ~1
  
  
  n.ind   <- ncol(trait)
  n.gene  <- nrow(trait)
  t0      <- t0.mki #apply(trait,1,function(x) t.test(x[gt1==1],x[gt1==0])$statistic)
  

  gt.perturb.mat <- NULL
  for (i in indSample){
    gt.temp     <-  gt1
    if(!is.na(gt1[i])){       #when NA, do not perturb     
      gt.temp[i]  <-  abs(gt1[i]-1) #  gt=1or2   gt1[i]-sign(gt1[i]-1.5)*1       
    }    
    gt.perturb.mat <- rbind(gt.perturb.mat, gt.temp)    #ith row is a genotype at the marker under studied that with jth RIL has been perturbed
  }
  dimnames(gt.perturb.mat) <- list(paste("gt1MK_",1:length(indSample),"_sample_perturbed",sep=""), colnames(trait))
  
  delta.t <- NULL      #n.ind x n.gene, flip each samples'gt(row), deltaT for each gene (column)at this particular mk    
  t.perturb <- tMatFunction (trait,gt.perturb.mat,fileName) #probe x nSample, here each colum corresponds t value when the gt of this sample has been perturbed
  delta.t <- abs(t.perturb) - matrix(abs(t0))[,rep(1,length(indSample))]
        
  #mean.delta.t    <- apply(delta.t,1,mean)
 rownames(delta.t)<- rownames(trait)
 colnames(delta.t) <- colnames(trait)[indSample]
  
  return(delta.t)     #delta.t mat
}



.deltaTMatAllmk_list <- function(gt,pheno,t.mat0=NULL,indSample=NULL,t.thres=3,fileName="test"){
#usuage: in this function, 
#gt can be full genotype data, then t.mat0 can be given or not(computed in this fucntion)
#if gt is only for permuted samples (used in permutation function, then t.mat0 must be input here which was computed based on all genotype data
#ind.sample decide for which sample a wls score is to be computed

    if(is.null(indSample)){
       indSample <- 1:ncol(pheno)                           #calculate the original deltaTAllmk_list
    }
    gt.samp        <- gt    
    pheno.samp     <- pheno
    
    if(is.null(t.mat0)){
      #t.mat0  <- t.mat.func(pheno.samp,gt.samp)
      t.mat0  <-tMatFunction (pheno.samp,gt.samp,fileName)
    }
    delta.t.mat.allmk   <- NULL
    for (i.mk in 1:nrow(gt.samp)){ #indSample, correct here July9 2012
        ind.sig0  <-  which(abs(t.mat0[,i.mk]) >  t.thres)
	
	      topN                     <- 25 #50
	      delta.t.1mk50            <- matrix(NA, topN ,length(indSample))   #correct here
        dimnames(delta.t.1mk50 ) <- list( paste("pheno",1:topN ,sep=""), colnames(gt.samp)[indSample] )   #a matrix: 50 sigGerne x RIL
        if(length(ind.sig0)   >1){
          if(length(ind.sig0)   >topN) {
            ind.sig <-  order(abs(t.mat0[,i.mk]),decreasing=T)[1:topN ] #sample(ind.sig0,50,replace=F)                     #now i use the top sig ones
          }
          else{ ind.sig <- ind.sig0}        
          pheno.sig               <- pheno.samp[ind.sig, ]    #correct here, no column index: indSample; correct2: pheno.samp instead of pheno
          gt1                     <- gt.samp[i.mk, ]          #correct here, no column index: indSample; correct2: gt.samp instead of gt
          t0.mki                  <- t.mat0[ind.sig,i.mk]
          #delta.t.1mk             <- deltaTmatAtOneMK(pheno.sig,gt1,t0.mki)   #flip each sample(colum), deltaT for each gene(row)             
          delta.t.1mk             <- .deltaTmatAtOneMK(pheno.sig,gt1,t0.mki,fileName=fileName,indSample=indSample)   #flip each sample(colum), deltaT for each gene(row)             
          delta.t.1mk50[1:length(ind.sig),] <- delta.t.1mk
	        rownames(delta.t.1mk50)[1:length(ind.sig)] <- rownames(pheno.sig)
         }  else   {                #no sig QTL at this mk, then no wls can be detected           
           delta.t.1mk50            <- delta.t.1mk50 
        } 
         
        delta.t.mat.allmk       <- c(delta.t.mat.allmk,list(delta.t.1mk50)) #length=nMK, the output element is always 50genes x length(sampled RILs)
    }
       
    return(delta.t.mat.allmk)    #a list with length=nMK, each element is a matrix nSample by nSigGene , 
# each element matrix is deltaT for each sigGene(row) when  sample i (column) is flipped
}


.area.my <- function(deltaT.1sampleP0){
  if(any(is.na(deltaT.1sampleP0))){
  		deltaT.1sampleP <- deltaT.1sampleP0[-which(is.na(deltaT.1sampleP0))]
 	} else{
  		deltaT.1sampleP <- deltaT.1sampleP0
 	}
	#this function calculate the area of t>0 in the psiv distribution=deltaT.1sampleP=deatlT value when one sample is pertubed(taking the complementary gt)
      #library(zoo)
      Y <-density(deltaT.1sampleP)
      #Y is density
      Avg.pos <- 0
      # construct lengths and heights            #http://stackoverflow.com/questions/3876219/calculating-an-area-under-a-continuous-density-plot
      xt <- diff(Y$x[Y$x>Avg.pos])
      yt <- rollmean(Y$y[Y$x>Avg.pos],2)
      # This gives you the area
      a <- sum(xt*yt)
      return(a)
}

print.wls <- function(x,...){
  .check.wlsObject(x)
  cat( "\n")
  cat("OUTPUT of reGenotyper is a list with several elements: \n\n")
  cat("[1] WLS score for each sample \n\n")
  print(x$wls.score)
  cat("\n[2] detected WLS sample names \n\n")
  print(x$wls.names)
  cat("\n[3] recovered optimal genotype \n\n")
  print(x$gt.opt)
  cat("\n[4] WLS p value for each sample from permutation \n\n")
  print(x$wls.pValue)
  cat("\n[5] WLS score for permuted sample \n\n " )  
  print(x$wls.score.permu)
  cat("\n[6] threshold used probability threshold to decide if a sample is mislabled based on permutation result:" )  
  print(x$wls.score.permu)
  cat( "\n")
}

plot.wls <- function(x,...){
  .check.wlsObject(x)
  par(mfrow=c(2,1))
  wls.ind       <- which(names(x$wls.score)%in%x$wls.names)
  mycol         <- rep("black",length(x$wls.score))
  mycol[wls.ind]<- "red"

  #first plot - wls score
  plot(x$wls.score, ylim=c(min(x$wls.score)*0.95, max(x$wls.score)*1.15),ylab="WLS score",xlab="Sample",pch=19,cex.lab=1.2,cex.axis=1.2,col=mycol,main="")
  if(length(wls.ind>0)){
    text(x=wls.ind,y=x$wls.score[wls.ind]*1.05,labels=x$wls.names,col="red",cex=1.1,pos=rep(c(1,3),ceiling(length(wls.ind)/2) ))
  }

  #second plot - wls probability
  plot(x$wls.pValue, ylim=c(0,1.15),ylab="Probability of being WLS",xlab="Sample",pch=19,cex.lab=1.2,cex.axis=1.2,col=mycol,main="")
  abline(h=x$thres,col="gray",lty=2)
  if(length(wls.ind>0)){
    text(x=wls.ind,y=x$wls.pValue[wls.ind]*1.05,labels=x$wls.names,col="red",cex=1.1,pos=rep(c(1,3),ceiling(length(wls.ind)/2) ))
  }
}

.check.wlsObject <- function(wlsObject){
  if(!any(class(wlsObject)=="wls")) stop("This is not an object of wls class\n")
  if(any(is.na(wlsObject$wls.score))) stop("This is not a correct object of wls class, it is missing $wls.score\n")
  if(any(is.na(wlsObject$wls.names))) stop("This is not a correct object of wls class, it is missing $wls.names\n")
  if(any(is.na(wlsObject$wls.pValue))) stop("This is not a correct object of wls class, it is missing $wls.pValue\n")
  if(any(is.na(wlsObject$wls.score.permu))) stop("This is not a correct object of wls class, it is missing $wls.score.permu\n")
  if(any(is.na(wlsObject$gt.opt))) stop("This is not a correct object of wls class, it is missing $gt.opt\n")
  if(any(is.na(wlsObject$thres))) stop("This is not a correct object of wls class, it is missing $thr\n")
}