################################################################################
# fstWC: a function co calculate weir and cockerhams fis, fit, and fst
################################################################################
fstWC<-function(x){
  badData <- sapply(x$indtyp, function(y){
    is.element(0, y)
  })
  if(sum(badData) > 0){
    nl <- x$nloci - (sum(badData))
  } else{
    nl <- x$nloci
  }
  gdData<-which(!badData)
  badData<-which(badData)
  if (nl == 1) {
    all_genot<-x$pop_list[[1]][,gdData]
    if(x$npops > 1){
      for(i in 2:x$npops){
        all_genot <- c(all_genot, x$pop_list[[i]][,gdData])
      }
    }
    all_genot <- matrix(all_genot, ncol = 1)
  } else {
    all_genot<-matrix(x$pop_list[[1]][,gdData], ncol = length(gdData))
    if(x$npops > 1){
      for(i in 2:x$npops){
        all_genot<-rbind(all_genot, x$pop_list[[i]][,gdData])
      }
    }
  }
  genot<-apply(all_genot,2,unique)
  genot<-lapply(genot, function(x){
    if (sum(is.na(x))>0){
      y<-which(is.na(x)==TRUE)
      x_new<-x[-y]
      return(x_new)
    } else {
      return(x)
    }
  })
  #count genotypes
  
  genoCount<-list()
  for(i in 1:ncol(all_genot)){
    genoCount[[i]]<-matrix(0,ncol=length(genot[[i]]))
    for(j in 1:length(genot[[i]])){
      genoCount[[i]][,j]<-length(which(all_genot[,i] == genot[[i]][j]))
    }
    if (x$gp==3){
      colnames(genoCount[[i]])<-paste(substr(genot[[i]],1,3),"/",
                                      substr(genot[[i]],4,6),sep="")
    } else if (x$gp==2){
      colnames(genoCount[[i]])<-paste(substr(genot[[i]],1,2),"/",
                                      substr(genot[[i]],3,4),sep="")
    }
  }
  
  h_sum<-list()
  for(i in 1:ncol(all_genot)){
    h_sum[[i]]<-vector()
    cnSplit<-strsplit(colnames(genoCount[[i]]),"/")
    for(j in 1:length(x$all_alleles[[gdData[i]]])){
      het_id1<-lapply(cnSplit, is.element, x$all_alleles[[gdData[i]]][j])
      het_id2<-lapply(het_id1, sum)
      het_id2<-as.vector(het_id2)
      het_id3<-which(het_id2==1)
      h_sum[[i]][j]<-sum(genoCount[[i]][1,het_id3])
    }
  }
  indtyp_tot<-lapply(x$indtyp, sum)
  kk_hsum <- lapply(1:ncol(all_genot), function(i){
    list(h_sum[[i]], indtyp_tot[[gdData[i]]])
  })
  kk_hbar<-lapply(kk_hsum, function(x){
    return(x[[1]]/x[[2]])
  })
  
  pdat <- lapply(1:ncol(all_genot), function(i){
    list(x$allele_freq[[gdData[i]]], x$indtyp[[gdData[i]]])
  })
  
  kk_p<-lapply(pdat, function(x){
    if(is.null(x[[1]])==FALSE){
      apply(x[[1]], 1, function(y){
        y*(2*x[[2]])
      })
    }
  })
  res<-matrix(0,(x$nloci+1),3)
  colnames(res)<-c("Fis_WC","Fst_WC","Fit_WC")
  rownames(res)<-c(x$loci_names, "All")
  A<-vector()
  a<-vector()
  b<-vector()
  c<-vector()
  for(i in 1:ncol(all_genot)){
    kknbar<-indtyp_tot[[gdData[i]]]/x$npops
    kknC<-(indtyp_tot[[gdData[i]]]-sum(x$indtyp[[gdData[i]]]^2)/
             indtyp_tot[[gdData[i]]])/(x$npops-1)
    kkptild<-kk_p[[i]]/(2*x$indtyp[[gdData[i]]])
    kkptild[kkptild=="NaN"]<-NA
    kkpbar<-colSums(kk_p[[i]])/(2*indtyp_tot[[gdData[i]]])
    kks2<-colSums(x$indtyp[[gdData[i]]]*
                    (kkptild-rep(kkpbar,each = x$npops))^2)/((x$npops-1)*kknbar)
    kkA<-kkpbar*(1-kkpbar)-(x$npops-1)*kks2/x$npops
    kka<-kknbar*(kks2-(kkA-(kk_hbar[[i]]/4))/(kknbar-1))/kknC
    kkb<-kknbar*(kkA-(2*(kknbar-1))*kk_hbar[[i]]/(4*kknbar))/(kknbar-1)
    kkc<-kk_hbar[[i]]/2
    A[i]<-sum(kkA)
    a[i]<-sum(kka)
    b[i]<-sum(kkb)
    c[i]<-sum(kkc)
    res[gdData[i],"Fis_WC"]<- round(1-sum(kkc)/sum(kkb+kkc),4)
    res[gdData[i],"Fst_WC"]<- round(sum(kka)/sum(kka+kkb+kkc),4)
    res[gdData[i],"Fit_WC"]<- round(1-sum(kkc)/sum(kka+kkb+kkc),4)
  }
  res[res=="NaN"]<-NA
  res[res==0.000]<-NA
  sumA<-sum(na.omit(A))
  suma<-sum(na.omit(a))
  sumb<-sum(na.omit(b))
  sumc<-sum(na.omit(c))
  res[(x$nloci+1),"Fis_WC"]<-round(1-sumc/(sumb+sumc),4)
  res[(x$nloci+1),"Fst_WC"]<-round(suma/(suma+sumb+sumc),4)
  res[(x$nloci+1),"Fit_WC"]<-round(1-sumc/(suma+sumb+sumc),4)
  #res[is.na(res)]<-NaN
  list(Fstats=res,
       multiLoc<-res[(x$nloci+1),])
}
################################################################################
# end fstWC
################################################################################