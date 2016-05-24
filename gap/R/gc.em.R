gc.em <- function(data, locus.label=NA, converge.eps=0.000001, maxiter=500, 
         handle.miss=0, miss.val=0, control=gc.control())
{
  if (control$xdata) {
     sex <- data[,1]
     data <- data[,-1]
  }
  for(p in c("haplo.stats")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!require(p, quietly = TRUE, character.only=TRUE))
        warning(paste("gc.em needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  tmp0<-geno.recode(data,miss.val=miss.val)
  geno<-tmp0$grec
  geno[is.na(geno)]<-0
  data<-as.matrix(geno)
  weight<-rep(1,dim(data)[1])
  nloci<-dim(data)[2]/2
  loci<-rep(0,nloci)
  for (i in 1:nloci)
  {
      loci[i]<-length(tmp0$alist[[i]]$allele) # max(data[,c(2*i-1,2*i)],na.rm=TRUE)
  }
  if(all(is.na(locus.label))) {
     locus.label<- paste("loc-",1:nloci,sep="")
  }
# to run genecounting
  if (control$xdata) data <- cbind(sex,data)
  data.gc<-genecounting(data,weight=weight,loci=loci,
           control=gc.control(xdata=control$xdata,
           eps=converge.eps,pl=0.001,maxit=maxiter,
           handle.miss=handle.miss,assignment=control$assignment,verbose=F))
  hap.prob<-data.gc$h
  hap.prob.noLD<-data.gc$h0
  lnlike<-data.gc$l1
  lr<-2*(data.gc$l1-data.gc$l0)
  df<-data.gc$npdat-sum(loci)-length(loci)
  niter<-data.gc$iter
  converge<-data.gc$converge
# to further extract information and obtain unique haplotypes
  hapas<-read.table(control$assignment)
  unlink(control$assignment)
  newnames<-c("subj","chr",locus.label,"post","hapid")
  names(hapas)<-newnames
  ncol<-nloci+4
  nrow<-dim(hapas)[1]/2
  indx1<-2*1:nrow-1
  indx2<-2*1:nrow
  indx.subj<-hapas$subj[indx1]
  hapdat<-hapas[,-c(1,2,ncol-1)]
  post<-hapas$post[indx1]
  hapid<-hapas$hapid
  one<-rep(1,nrow*2)
  hapdat<-cbind(hapdat,one)
  with(hapdat,{
  tmp<-by(hapdat,one,unique)
  haplotype<-as.matrix(tmp[[1]])
  tmp<-order(haplotype[,nloci+1])
  haplotype<-haplotype[tmp,1:(dim(haplotype)[2]-2)]
  dimnames(haplotype)<-list(1:length(haplotype[,1]),locus.label)
  hap1<-hapid[indx1]
  hap2<-hapid[indx2]
  uhap <- sort(unique(hapid))
  hap.prob<-hap.prob[uhap]
  nreps<-tapply(indx.subj,indx.subj,length)
# 13/11/2003
# haplotype trend regression
# assign.dat already has sequential number to avoid duplicate IDs
  idx.subj<-sort(unique(indx.subj))
  N<-length(idx.subj)
  P<-length(uhap)
  idx.subj<-cbind(1:N,idx.subj)
  idx.uhap<-cbind(1:P,uhap)
  htrtable<-matrix(rep(0,N*P),nrow=N)
  for(l in 1:nrow)
  {
    i<-idx.subj[,1][idx.subj[,2]==indx.subj[l]]
    j1<-idx.uhap[,1][idx.uhap[,2]==hap1[l]]
    htrtable[i,j1]<-htrtable[i,j1]+post[l]
    j2<-idx.uhap[,1][idx.uhap[,2]==hap2[l]]
    htrtable[i,j2]<-htrtable[i,j2]+post[l]
  }
  htrtable<-htrtable/2
  dimnames(htrtable)<-list(NULL,as.character(uhap))
  list(lnlike=lnlike,lr=lr,
       hap.prob=hap.prob,hap.prob.noLD=hap.prob.noLD,indx.subj=indx.subj,
       post=post,hap1code=hap1,hap2code=hap2,haplotype=grec2g(haplotype,nloci,tmp0),
       nreps=nreps,converge=converge,niter=niter,uhap=uhap,htrtable=htrtable)
  })
}
