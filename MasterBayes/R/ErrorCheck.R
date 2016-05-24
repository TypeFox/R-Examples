"ErrorCheck"<-function(sP, tP, pP, PdP, GdP, unique_id, nbeta, nlinked){

  if(is.null(PdP$data)){Ppresent<-FALSE}else{Ppresent<-TRUE}
  if(is.null(GdP$G)){Gpresent<-FALSE}else{Gpresent<-TRUE}

  if(length(unique(PdP$id))!=length(PdP$id)){ 
    if(length(unique(PdP$timevar))<2){
      stop("multiple records per individual but timevar does not vary")
    }
  }

  if(Ppresent==TRUE){
    if(length(PdP$USdam)==1 & PdP$USdam[1]==FALSE){
      sP$estUSdam<-FALSE
      nusd<-0
    }else{
      nusd<-length(unique(PdP$USdam))
    }

    if(length(PdP$USsire)==1 & PdP$USsire[1]==FALSE){
      sP$estUSsire<-FALSE
      nuss<-0
    }else{
      nuss<-length(unique(PdP$USsire))
    }

    if(sum(nbeta)==0){sP$estbeta<-FALSE}
  }

  No.E<-length(unique(GdP$categories))*length(GdP$G)*GdP$perlocus+((1-GdP$perlocus)*length(unique(GdP$categories)))

  if(Gpresent==FALSE){
    sP$estG<-FALSE
    sP$estE1<-FALSE
    sP$estE2<-FALSE
    sP$estA<-FALSE
    sP$estP<-FALSE
  }

  if(sP$estG==FALSE & is.null(sP$G)){
    sP$estE1<-FALSE
    sP$estE2<-FALSE
    sP$estA<-FALSE
  }

  if(Ppresent==FALSE){
    sP$estP<-FALSE
    sP$estbeta<-FALSE
    sP$estUSsire<-FALSE
    sP$estUSdam<-FALSE
  }

  if(is.null(PdP$sex)==FALSE){
    if(sum(names(table(PdP$sex))%in%c("Male", "Female")==FALSE)>0){
      stop("sex levels != Male Female")
    }
  }  
############### check starting parameteristations if they exist ###############################################

  if(is.null(sP$id)==FALSE & Ppresent==TRUE){
    if(FALSE%in%(sP$id%in%unique_id)){
      stop("indivdiuals exist in sP$id not in PdataPed object")
    }
    if(FALSE%in%(unique_id%in%sP$id)){
      stop("indivdiuals exist in PdataPed object but not in sP$id")
    }
  }

  if(is.null(sP$G)==FALSE){
    if(is.genotype(sP$G[[1]])==FALSE){
      sP$G<-genotype.list(sP$G)
      if(length(sP$id)!=length(sP$G[[1]][,1])){
        stop("different number of starting genotypes than in sP$id")
      }
    }
    if(is.null(GdP$G)==FALSE){
      if(length(GdP$G)!=length(sP$G)){
        stop("different number of loci in starting genotypes and genotype data")
      }
    }
  }
  if(is.null(sP$A)==FALSE){
    if(length(names(sP$A[[1]]))==0){
      stop("starting allele frequencies do not have allele names")
    }
    if(is.null(sP$G)==FALSE){
      if(length(sP$A)!=length(sP$G)){
        stop("different number of loci in starting genotypes and starting allele frequencies")
      }
      for(i in 1:length(sP$A)){
        if(any((allele.names(sP$G[[i]])%in%names(sP$A[[i]]))==FALSE)){
          stop("some alleles in sP$G are not in sP$A")
        }
      }
    }
    if(is.null(GdP$G)==FALSE){
      if(length(sP$A)!=length(GdP$G)){
        stop("different number of loci in starting allele frequencies and genotype data")
      }
      for(i in 1:length(sP$A)){
        if(any((allele.names(GdP$G[[i]])%in%names(sP$A[[i]]))==FALSE)){
          stop("some alleles in GdP$G are not in sP$A")
        }
      }
    }
  }

  if(is.null(sP$E1)==FALSE){
    if(length(sP$E1)==1){
      sP$E1<-rep(sP$E1, No.E)
    }
    if(length(sP$E1)!=No.E){
      stop("number of error categories does not match length of sP$E1")
    } 
  }
 
  if(is.null(sP$E2)==FALSE){
    if(length(sP$E2)==1){
      sP$E2<-rep(sP$E2, No.E)
    }
    if(length(sP$E2)!=No.E){
      stop("number of error categories does not match length of sP$E2")
    } 
  }

  if(is.null(sP$beta)==FALSE){
    if(length(sP$beta)!=sum(nbeta[1:6])){
      stop("number of parameters to estimate does not match length of sP$beta")
    } 
  }

  if(is.null(sP$USsire)==FALSE){
    if(length(sP$USsire)==1){
      sP$USsire<-rep(sP$USsire, nuss)
    }
    if(length(sP$USsire)!=nuss){
      stop("number of unsampled sire categories does not match length of sP$USsire")
    } 
  }

  if(is.null(sP$USdam)==FALSE){
    if(length(sP$USdam)==1){
      sP$USdam<-rep(sP$USdam, nusd)
    }
    if(length(sP$USdam)!=nusd){
      stop("number of unsampled dam categories does not match length of sP$USsire")
    } 
  }

  if(is.null(sP$ped)==FALSE){
    if(FALSE%in%(unique_id%in%sP$ped[,1])){
      if(FALSE%in%(unique_id[PdP$id[which(PdP$offspring==1)]]%in%sP$ped[,1])){
        stop("sP$ped must have rows for all offspring or all indiviuals")
      }else{
        base<-setdiff(unique_id, sP$ped[,1]) 
        sP$ped<-rbind(cbind(as.character(base), rep(NA,length(base)), rep(NA,length(base))),  sP$ped)      
      }
    }
    if(FALSE%in%(sP$ped[,1]%in%unique_id)){
      stop("some individuals in the id column of sP$ped are not in PdP")
    }
    if(FALSE%in%(na.omit(sP$ped[,2])%in%unique_id)){
      stop("some dams in sP$ped are not in PdP")
    }
    if(FALSE%in%(na.omit(sP$ped[,3])%in%unique_id)){
      stop("some sires in sP$ped are not in PdP")
    }
    if(is.null(PdP$sex)==FALSE){
      if("Male"%in%PdP$sex[match(match(na.omit(sP$ped[,2]), unique_id), PdP$id)]){
       stop("some starting dams are recorded as males")
      }
      if("Female"%in%PdP$sex[match(match(na.omit(sP$ped[,3]), unique_id), PdP$id)]){
       stop("some starting sires are recorded as females")
      }                 
    }
  }


###################################################################################################################
######################## checking whether things to be estimated can be ########################################### ###################################################################################################################

  if(sP$estP==TRUE){
    if(sP$estG==FALSE & is.null(sP$G)){
      print("using an approximation for genotyping error")
    }
  }

  if(Gpresent==TRUE & Ppresent==TRUE){
    if(FALSE%in%(GdP$id%in%PdP$id)){
      stop("genotype data exists for individuals not in PdataPed object")
    }
    if(FALSE%in%(PdP$id%in%GdP$id)){
      stop("some individuals in PdataPed object have no genotype data: replace with NA")
    }
  }

  if(sP$estbeta==TRUE){
    if(sP$estP==FALSE & (sum(nbeta[1:2])>0 | sum(nbeta[5:6])>0) & is.null(sP$ped)){
      stop("beta cannot be estimated because dams are not given or cannot be estimated")
    }
    if(sP$estP==FALSE & (sum(nbeta[3:4])>0 | sum(nbeta[5:6])>0) & is.null(sP$ped)){
      stop("beta cannot be estimated because sire are not given or cannot be estimated")
    }
  }

  if(sP$estP==TRUE & sP$estA==FALSE & sP$estG==FALSE & is.null(sP$A)){
      warning("allele freqencies estimated from GdP$G")
    }

  if(sP$estE1==TRUE){
    if(sP$estP==FALSE & is.null(sP$ped) & length(GdP$id)==length(unique(GdP$id)) & pP$E1[1]==999){
      stop("E1 cannot be estimated because pedigree cannot be estimated and duplicate samples do not exist")
    }
  }

  if(sP$estE2==TRUE){
    if(sP$estP==FALSE & is.null(sP$ped) & length(GdP$id)==length(unique(GdP$id)) & pP$E1[1]==999){
      stop("E2 cannot be estimated because pedigree cannot be estimated and duplicate samples do not exist")
    }
  }

  if(sP$estA==TRUE){
    if(is.null(GdP$G) & is.null(sP$G)){
      stop("alelle frequnecies cannot be estimated because genetic data do not exist")
    }
  }

########## checking whether prior specifications are valid ##########################################################

 if(pP$E1[1]!=999){
   if(dim(pP$E1)[1]!=No.E){
     stop("Matrix specifying prior distribution of E1 is the wrong dimension")
   }
   if(dim(pP$E1)[2]!=2){
      stop("Matrix specifying prior distribution of E1 is the wrong dimension")     
   }
 }

 if(pP$E2[1]!=999){
   if(dim(pP$E2)[1]!=No.E){
     stop("Matrix specifying prior distribution of E2 is the wrong dimension")
   }
   if(dim(pP$E2)[2]!=2){
      stop("Matrix specifying prior distribution of E2 is the wrong dimension")     
   }
 }

 if(sP$estUSsire==TRUE & sP$estUSdam==TRUE){
   if((pP$USsire$mu[1]!=999 & pP$USdam$mu[1]==999) | (pP$USsire$mu[1]==999 & pP$USdam$mu[1]!=999)){
     stop("unsampled mothers and fathers are updated as a block and prior speicifictaions must exist for both, or neither")
   }
 }
 if(pP$beta$mu[1]!=999 | pP$beta$sigma[1]!=999){
   if(length(pP$beta$mu)!=(sum(nbeta)-nlinked)){
     stop("Vector specifying the mean of the prior distribution for beta is the wrong length")
   }
   if(dim(pP$beta$sigma)[1]!=(sum(nbeta)-nlinked)| dim(pP$beta$sigma)[2]!=(sum(nbeta)-nlinked)){
     stop("Matrix specifying the (co)variance of the prior distribution for beta is the wrong dimension")
   }
 }

 if(pP$USdam$mu[1]!=999 | pP$USdam$sigma[1]!=999){
   if(length(pP$USdam$mu)!=nusd){
     stop("Vector specifying the mean of the prior distribution for USdam is the wrong length")
   }
   if(length(pP$USdam$sigma)!=nusd){
     stop("Vector specifying the variance of the prior distribution for USdam is the wrong length")
   }
 }

 if(pP$USsire$mu[1]!=999 | pP$USsire$sigma[1]!=999){
   if(length(pP$USsire$mu)!=nuss){
     stop("Vector specifying the mean of the prior distribution for USsire is the wrong length")
   }
   if(length(pP$USsire$sigma)!=nuss){
     stop("Vector specifying the variance of the prior distribution for USsire is the wrong length")
   }
 }
sP
}
 
