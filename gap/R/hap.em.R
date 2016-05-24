hap.em <- function(id,data,locus.label=NA,converge.eps=0.000001,maxiter=500,miss.val=0)
{
  for(p in c("haplo.stats")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!require(p, quietly = TRUE, character.only=TRUE))  
        warning(paste("hap.em needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  tmp0 <- geno.recode(data,miss.val=miss.val)
  geno <- as.matrix(tmp0$grec)
  geno[is.na(geno)] <- 0
  data<-as.matrix(geno)
  nloci<-dim(data)[2]/2
  loci<-rep(0,nloci)
  for (i in 1:nloci)
  {
     loci[i]<-length(tmp0$alist[[i]]$allele) #  max(data[,c(2*i-1,2*i)],na.rm=TRUE)
  }
  if(all(is.na(locus.label))) {
     locus.label<- paste("loc-",1:nloci,sep="")
  }
##
  zz <- hap.control(ss=1)
  z<- hap(id=id,data=data,nloci,loci=loci,control=zz)
#
  tmp1<-read.table(zz$hapfile,header=T)
# unlink("hap.out")
  haplotype<-as.matrix(tmp1[,1:nloci])
  dimnames(haplotype)<-list(1:length(haplotype[,1]),locus.label)
  uhap<-hapid<-1:(dim(tmp1)[1])
  tmp1<-data.frame(tmp1,hapid)
  hap.prob<-tmp1[,nloci+1]
#
  tmp2<-read.table(zz$assignfile,header=T)
# unlink("assign.out")
  nrow<-dim(tmp2)[1]/2
  indx1<-2*1:nrow-1
  indx2<-2*1:nrow
  indx.subj<-tmp2[indx1,1]
  post<-tmp2[indx1,nloci+3]

  tmp<-merge(tmp1[,-(nloci+1)],tmp2[,-c(1,2)],sort=F)
  hap1<-tmp[indx1,(nloci+1)]
  hap2<-tmp[indx2,(nloci+1)]
  nreps<-tapply(indx.subj,indx.subj,length)

  list (lnlike=z$l1,hap.prob=hap.prob,indx.subj=indx.subj,post=post,
        hap1code=hap1,hap2code=hap2,haplotype=grec2g(haplotype,nloci,tmp0),nreps=nreps,
        converge=z$converge,niter=z$niter,uhap=uhap)
}

# 15-10-03 WAH
# 16-10-03 UCL office to add sort=F in merge
# 25-09-04 first attempt to robust label-handling
# 29-07-07 add library(haplo.stats)
