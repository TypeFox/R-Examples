# $Id: hapmcmc.R 1352 2012-08-14 14:21:35Z warnes $
#
# Code contributed by David Duffy <davidD@qumr.edu.au>:
#
# "If you are interested, this is a toy/prototype for haplotyping via MCMC.  
#  It is much slower than Dan Schaid's haplo.em, but does give the same
#  answers ;)"

#
# Routines for handling genotypes
#
# Convert "1/2" to 1,2
#
geno.as.array <- function(genotypes,renumber=FALSE,miss=NULL,gtp.sep="/") {
  mknum<-function(genotypes, renumber=FALSE, gtp.sep="/") {
    alleles<- strsplit(genotypes, gtp.sep)
    gtp<-cbind(sapply(alleles, function(x) x[1], simplify=TRUE),
               sapply(alleles, function(x) x[2], simplify=TRUE))
    if (renumber) {
      alleles<-unique(unlist(alleles))
      gtp[,1]<-as.numeric(factor(gtp[,1],levels=alleles))
      gtp[,2]<-as.numeric(factor(gtp[,2],levels=alleles))
    }
    if (is.null(miss)) {
      gtp[is.na(genotypes),]<-NA
    }else{
      gtp[is.na(genotypes),]<-miss
    } 
    gtp
  } 
  if (is.null(ncol(genotypes)) || ncol(genotypes)==1) {
    res<-mknum(genotypes, renumber=renumber)
  }else{
    res<-data.frame(mknum(genotypes[,1], renumber=renumber))
    for(i in 2:ncol(genotypes)) {
      res<-cbind(res,mknum(genotypes[,i], renumber=renumber))
    } 
    colnames(res)<-c(t(outer(names(genotypes),1:2,paste,sep="."))) 
  } 
  apply(res,2,as.character)
} 
#
hap <- function(genotypes) {
  res<-geno.as.array(genotypes)
  nc<-ncol(res)
  hap1<-res[,seq(1,nc,2)]
  hap2<-res[,seq(2,nc,2)]
  loci<-colnames(genotypes)
  colnames(hap2)<-colnames(hap1)<-loci
  list(hap1=hap1, hap2=hap2, class="haplotype")
}
hapshuffle <- function(haplotypes, hfreq=NULL, ambiguous=NULL, verbose=FALSE, set) {
  if (is.null(hfreq)) hfreq<-hapfreq(haplotypes, set=set)  
  if (is.null(ambiguous)) ambiguous<-hapambig(haplotypes)
  nloci<-ncol(haplotypes$hap1)
  nobs<-nrow(haplotypes$hap1)
  for(ind in ambiguous) {
    prop<-curr<-list(hap1=haplotypes$hap1[ind,], hap2=haplotypes$hap2[ind,])
    swap<-sample(c(TRUE,FALSE),nloci,replace=TRUE)
    if (any(swap)) {
      tmp<-prop$hap1[swap]
      prop$hap1[swap]<-prop$hap2[swap]
      prop$hap2[swap]<-tmp
    }
    o1<-paste(curr$hap1,collapse=":")  
    o2<-paste(curr$hap2,collapse=":")  
    n1<-paste(prop$hap1,collapse=":")  
    n2<-paste(prop$hap2,collapse=":")  
    pos.o1<-match(o1,names(hfreq))
    pos.o2<-match(o2,names(hfreq))
    pos.n1<-match(n1,names(hfreq))
    pos.n2<-match(n2,names(hfreq))
    pn<-(hfreq[pos.n1]+0.5)*(hfreq[pos.n2]+0.5)
    po<-(hfreq[pos.o1]+0.5)*(hfreq[pos.o2]+0.5)
    qa<-pn/po
    if (verbose) cat("Person ",ind," ",qa," ",o1,"/",o2," -> ",n1,"/",n2,sep="")
    if (qa>runif(1)) {
      if (verbose) cat(" Accepted\n")
      haplotypes$hap1[ind,]<-prop$hap1
      haplotypes$hap2[ind,]<-prop$hap2
      hfreq[pos.n1]<-hfreq[pos.n1]+1
      hfreq[pos.n2]<-hfreq[pos.n2]+1
      hfreq[pos.o1]<-hfreq[pos.o1]-1
      hfreq[pos.o2]<-hfreq[pos.o2]-1
    }else if (verbose) {
      cat(" Unchanged\n")
    }
  }
  list(hfreq=hfreq, haplotypes=haplotypes, class="hapshuffle")
}

hapambig <- function(haplotypes) {
  which(apply(haplotypes$hap1!=haplotypes$hap2,1,sum)>1)
}

hapenum <- function(haplotypes) {
  dat<-rbind(haplotypes$hap1, haplotypes$hap2)
  dat<-dat[complete.cases(dat),]
  set<-unique(dat[,1])
  for(i in 2:ncol(dat)) set<-outer(set,unique(dat[,i]),paste,sep=":") 
  factor(set)
}

hapfreq <- function(haplotypes, set=NULL) {
  if (is.null(set)) set<-hapenum(haplotypes)
  hap1<-apply(haplotypes$hap1[complete.cases(haplotypes$hap1),],1,paste,collapse=":")
  hap2<-apply(haplotypes$hap2[complete.cases(haplotypes$hap2),],1,paste,collapse=":")
  dat<-c(hap1,hap2)
  table(factor(dat,levels=set))
}

hapmcmc <- function(gtp, B=1000) {
  tot<-2*nrow(gtp)
  hap.dat<-hap(gtp)
  hap.set<-hapenum(hap.dat)
  hap.amb<-hapambig(hap.dat)
  hap.new<-list(hfreq=hapfreq(hap.dat, set=hap.set), haplotypes=hap.dat)
  res<-matrix(nrow=B, ncol=length(hap.set))
  colnames(res)<-as.character(hap.set)
  rownames(res)<-1:B
  for(i in 1:B) {
    hap.new<-hapshuffle(hap.new$haplotypes,hfreq=hap.new$hfreq,ambiguous=hap.amb, set=hap.new$set)
    res[i,]<-hap.new$hfreq
  }
  apply(res,2,mean)/tot
}

mourant <- function(n) {
  tab<-matrix(c(91,32,5,147,78,17,85,75,7), nrow=3)
  rownames(tab)<-c("M/M","M/N","N/N")
  colnames(tab)<-c("S/S","S/s","s/s")
  dat<-as.data.frame.table(tab)
  p<-dat$Freq/sum(dat$Freq)
  dat[sample(1:nrow(dat),n,replace=TRUE,prob=p),1:2]
}
