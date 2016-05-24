sim.control<-function(
         simulate=FALSE, # do we want to simulate?
         # relationships to simulate -- by default, don't do cousins
         rships=c("unrelated","MZtwins","parent-offspring","full-sibs","half-sibs"),
         nsim=rep(200,length(rships)), # number of gene drops per rship
         userdat=NULL,
         geno.err=1/1000,
         hom2hom.err=0,
         fitLD=TRUE, #fit LD models?
         LDfiles=NULL, # list of LD files
         cl=NULL) #SNOW cluster object
{
  names(nsim)<-rships
  if(!simulate) {
    rships<-nsim<-geno.err<-hom2hom.err<-fitLD<-LDfiles<-cl<-NULL
  }
  list(simulate=simulate,rships=rships,nsim=nsim,userdat=userdat,
       geno.err=geno.err, hom2hom.err=hom2hom.err,fitLD=fitLD, 
       LDfiles=LDfiles,cl=cl)
}

################################################################
simIBD=function(snpmatlist,simparams){

    # Lots of decisions to be made in this function:
    # - User wants LD?
    #   * If yes, do we have to fit LD models?
    # But first, since both LD fitting and gene drops are done separately
    # for each chromosome. Create a list of data needed for each chrom.

    dat.by.chrom=split.by.chrom(snpmatlist,simparams)
    if(simparams$fitLD) {
      if(is.null(simparams$LDfiles)) {
        # Need to fit the LD models and save in parfiles
        GMLDfct(dat.by.chrom,simparams)
      } # else do nothing, because LD model files already exist
    } else {
      # User doesn't want LD model -- write regular parfiles
      lapply(dat.by.chrom,parfile.noLD)
    }
    snpmats<-genedrops(dat.by.chrom,simparams)
    return(snpmats)
}

split.by.chrom=function(snpmatlist,simparams){
        
    popsam=snpmatlist$subject.support$popsam
    snpobjr=snpmatlist$snp.data[popsam,]
    genloc=snpmatlist$snp.support$Gen_loc
    chrvec=snpmatlist$snp.support$Chromosome
        
    # We used to set up two lists: snpobjlr and genlocl.
    # Changed recently to set up one list called gdat with all info in it.
    # This makes it easier to re-write the code to use a cluster.
    gdat<-list(); c.count<-0
    for (i in 1:22){
        if(is.element(i,chrvec)=="TRUE"){
            c.count<-c.count+1 #found another chromosome
            rslist=snpobjr[,chrvec==i]
            loclist=genloc[chrvec==i]
            if(simparams$fitLD) { #parfile for gene drops will include LD model
              if(!is.null(simparams$LDfiles)) {
                GDparfile=simparams$LDfiles[c.count]
              } else {
                GDparfile=paste("GD",i,".ld.par",sep="")
              }
            } else {
              GDparfile=paste("GD",i,".par",sep="")
            }
            gdat[[c.count]]<-list(snp.data=rslist,genmap=loclist,chrom=i,
                                      GDparfile=GDparfile)
        }
    }
        
    nl=length(gdat)
    #ordering the snps according to genetic map position
    for (j in 1:nl){
        snpchr=gdat[[j]]$snp.data
        locsnps=as.numeric(as.vector(gdat[[j]]$genmap))
        ordsnps=order(locsnps)
        gdat[[j]]$snp.data=snpchr[,ordsnps]
        gdat[[j]]$genmap=locsnps[ordsnps]
    }
    return(gdat)
}

GMLDfct=function(gdat,simparams){

     if(is.null(simparams$cl)) { # no cluster: use lapply 
       unlist(lapply(gdat,parfile.fitLD))
     } else { # can use clusterApplyLB to distribute work over cluster
       require(parallel)
       unlist(clusterApplyLB(simparams$cl,gdat,parfile.fitLD))
     }
}
  
parfile.fitLD<-function(gdat) {
  # Use Thomas' FitGMLD to fit an LD model for founder haplotypes.
  # Computations become too much  of a burden for more than about 100 
  # subjects, so subset if necessary.
  if(nrow(gdat$snp.data)>100) {
    s.ind<-sample(1:nrow(gdat$snp.data))
    gdat$snp.data<-gdat$snp.data[s.ind,]
  }
  fp<-paste("flipped",gdat$chrom,".ped",sep="")
  write.pedfile(pedinf="unrelated",snp.data=gdat$snp.data,
        file=fp,transpose=TRUE)
  pf<-paste("parfile",gdat$chrom,".par",sep="")
  write.parfile(gdat$snp.data,gdat$genmap,file=pf)
  FitGMLD(pf,fp,gdat$GDparfile)
  unlink(c(pf,fp))
  return(gdat$GDparfile) #Will need to output names of parfiles with LD models
}     
  
parfile.noLD<-function(gdat) {
  # No fitted LD model, so just write a regular LINKAGE parfile.
  write.parfile(gdat$snp.data,gdat$genmap,file=gdat$GDparfile)
}
    
#######################
### gendrop simulations

genedrops<-function(gdat,simparams){

    # Do gene drops for all relationships of interest
    ursnpmat<-mzsnpmat<-posnpmat<-fssnpmat<-hssnpmat<-cosnpmat<-usersnpmat<-NULL
    if("unrelated" %in% simparams$rships) {
      ursnpmat=gendrp(gdat,peddat=ur.peddat(),
                  simparams$nsim["unrelated"],simparams)
    }
    if("MZtwins" %in% simparams$rships) {
     # Won't do gene drops, just sample subjects and apply genotyping error
      mzsnpmat=mztwinsim(gdat, simparams$nsim["MZtwins"],simparams)
    }
    if("parent-offspring" %in% simparams$rships) {
      posnpmat=gendrp(gdat,peddat=po.peddat(),
                  simparams$nsim["parent-offspring"],simparams)
    }
    if("full-sibs" %in% simparams$rships) {
      fssnpmat=gendrp(gdat,peddat=fs.peddat(),
                  simparams$nsim["full-sibs"],simparams)
    }
    if("half-sibs" %in% simparams$rships) {
      hssnpmat=gendrp(gdat,peddat=hs.peddat(),
                  simparams$nsim["half-sibs"],simparams)
    }
    if("cousins" %in% simparams$rships) {
      cosnpmat=gendrp(gdat,peddat=co.peddat(),
                  simparams$nsim["cousins"],simparams)
    }
    if("user" %in% simparams$rships) {
      usersnpmat=gendrp(gdat,peddat=user.peddat(simparams$userdat),
                  simparams$nsim["user"],simparams)
    }
    snpmats=list(mz=mzsnpmat,ur=ursnpmat,po=posnpmat,fs=fssnpmat,hs=hssnpmat,co=cosnpmat,user=usersnpmat,
                 LDfiles=findLDfiles(gdat))#Store LDfiles for output from simIBD
    return(snpmats)
}
    
mztwinsim=function(gdat,nreps,simparams){
  if(length(gdat) ==1) {
    snpmat<-doTwinsim(gdat[[1]],nreps,simparams)
  } else {
    snpmatlist<-lapply(gdat,doTwinsim,nreps,simparams)
    # cbind results into one big snp.matrix object
    snpmat<-Reduce(cbind,snpmatlist)
  }
  return(snpmat)
}


    
gendrp=function(gdat,peddat,nreps,simparams){

    if(is.null(simparams$cl)) {
      snpmatlist<-lapply(gdat,doGeneDrops,peddat,nreps,simparams)
    } else {
      require(parallel)
      snpmatlist<-clusterApply(simparams$cl,gdat,doGeneDrops,peddat,nreps,simparams)
    }
        
    # cbind results into one big snp.matrix object
    nl=length(snpmatlist)
    snpmat<-snpmatlist[[1]]
    if(nl > 1){
        for (j in 2:nl){
            snpmat=cbind(snpmat,snpmatlist[[j]])
        }
    }
        
    return(snpmat)
}
doTwinsim<-function(gdat,n,simparams) {
  # First draw a sample of n, with replacement, from the study data.
  s.inds <- sample(1:nrow(gdat$snp.data),n,replace=TRUE)
  # Replicate each sample to make it a pair.
  so1<-gdat$snp.data[s.inds,]
  row.names(so1)<-as.character(1:n)
  so2<-gdat$snp.data[s.inds,]
  row.names(so2)<-as.character((n+1):(2*n))
  # Convention for output is all pairmember1's followed by all pairmember2's.
  snpout<-rbind(so1,so2)
  # Apply genotyping errors
  nc<-ncol(gdat$snp.data)
  for(i in 1:nrow(snpout)) {
    # Sample genotypes to have errors with rate geno.err
    bb<-as.logical(rbinom(nc,size=1,prob=simparams$geno.err))
    snpout[i,bb]<-apply.errs(snpout[i,bb],simparams$hom2hom.err)
  }
  return(snpout)
}

# Worker function to doGeneDrops to do n gene drops given the pedigree 
# information. This function is called by the functions specific to each 
# relationship of interest
doGeneDrops<-function(gdat,pedinfo,n,simparams) {

  # We simulate complete data on the pedigree and then impose (i) genotyping
  # errors and (ii) missing
  # data patterns of randomly sampled pairs of subjects. GeneDrops needs
  # a LINKAGE pedfile as input. The pedigree information 
  # for a regular LINKAGE pedfile (and that's it) is in pedinfo. 
  # To complete the pedfile, we need to add marker information to this 
  # pedigree info.  As GeneDrops to simulate complete data on all subjects, 
  # the marker data we pass to GeneDrops is not actually used. Therefore, 
  # an empty (all missings) snp.matrix object of the right dimension will 
  # work as marker info.
  # The snp.matrix missing data code can be obtained by coercing 0
  # to the data-type "raw" with as.raw(). Use as.raw(0) to create the 
  # missing data code for snp.matrix objects. To create a snp.matrix object
  # of missing values, the command is:
  nr<-nrow(pedinfo); nc<-ncol(gdat$snp.data)
  sdat<-new("snp.matrix",data=as.raw(0),nrow=nr, ncol=nc,
             dimnames=list(as.character(1:nr),as.character(1:nc)))


  # Create input pedfile
  inped<-paste("in",gdat$chrom,".ped",sep="")
  write.pedfile(pedinfo,sdat,inped, sep=" ", eol="\n", na="0")
  # Do gene drops and save in an output pedfile
  outped<-paste("GDout",gdat$chrom,".ped",sep="")
  GeneDrops(gdat$GDparfile,inped,n,complete.data=TRUE, outped)
  # GeneDrops can fail but not report any error -- the output file may 
  # not exist.
  if(!file.exists(outped)) {
    warning("Gene drops failed for chromosome",gdat$chrom,"-- GeneDrops java program \n may have run out of memory and IBD coefficient clusters may not be valid.\n See the Note in the GeneDrops help file for instructions on how to increase java heap space.")
    return(NULL)
  }
  snpout=suppressWarnings(read.pedfile(file=outped,snp.names=colnames(gdat$snp.data))$snp.data)
  unlink(inped)
  unlink(outped)
  row.names(snpout)<-as.character(1:nrow(snpout))
  # Keep just the two subjects of interest, having ID's 1 or 2 (IDs are stored
  # in col 2 of pedinfo). Annick want's all ID=1's first, then ID=2's
  snpout<-rbind(snpout[pedinfo[,2]==1,],snpout[pedinfo[,2]==2,])
  #pairs<-rep(1:2,n)
  #snpout<-rbind(snpout[pairs==1,],snpout[pairs==2,])

  # Apply genotyping errors
  for(i in 1:nrow(snpout)) {
    # Sample genotypes to have errors with rate gerr
    bb<-as.logical(rbinom(nc,size=1,prob=simparams$geno.err))
    snpout[i,bb]<-apply.errs(snpout[i,bb],simparams$hom2hom.err)
  }

  # Select missing data patterns from study subjects and impose these on the 
  # simulated data.
  # First draw a random sample of n pairs from the study data.
  s.inds <- sample(1:nrow(gdat$snp.data),2*n,replace=TRUE)
  # Now loop over simulated data and impose missingness patterns from sampled
  # pairs
  for(i in 1:(2*n)) {
    # Need to test for the missing value in a snp.matrix object.
    miss.ind<-gdat$snp.data[s.inds[i],] == as.raw(0)
    snpout[i,miss.ind]<-as.raw(0)
  }

  return(snpout)
}

# Functions to set up pedigree information part of pedfiles
po.peddat<-function() { # Parent-offspring
    pedsize <- 3
    pedids <- rep(1,pedsize) 
    ids <- (1:pedsize) 
    dadids <- c(0,3,0)
    momids <- c(0,1,0)
    zeros <- rep(0,pedsize)
    genders <- c(2,2,1)
    ones <- rep(1,pedsize) #9th column(1 for the starting ind of the pedigree)
    peddat <- cbind(pedids,ids,dadids,momids,zeros,zeros,zeros,genders,ones)
    return(peddat)
}

fs.peddat=function(){ # Full sibs
    pedsize <- 4
    pedids <- rep(1,pedsize)
    ids <- (1:pedsize)
    dadids <- c(3,3,0,0)
    momids <- c(4,4,0,0)
    zeros <- rep(0,pedsize)
    genders <- c(2,2,1,2)
    ones <- rep(1,pedsize)
    peddat <- cbind(pedids,ids,dadids,momids,zeros,zeros,zeros,genders,ones)
    return(peddat)
}

hs.peddat=function(){ # Half sibs
        pedsize <- 5
        pedids <- rep(1,pedsize) 
        ids <- (1:pedsize) 
        dadids <- c(3,5,0,0,0)
        momids <- c(4,4,0,0,0)
        zeros <- rep(0,pedsize)
        genders <- c(2,2,1,2,1)
        ones <- rep(1,pedsize)
        peddat <- cbind(pedids,ids,dadids,momids,zeros,zeros,zeros,genders,ones)
        return(peddat)
}

ur.peddat=function(){ # unrelated
        pedsize <- 2
        pedids <- rep(1,pedsize) 
        ids <- (1:pedsize) 
        dadids <- c(0,0)
        momids <- c(0,0)
        zeros <- rep(0,pedsize)
        genders <- c(2,2)
        ones <- rep(1,pedsize)
        peddat <- cbind(pedids,ids,dadids,momids,zeros,zeros,zeros,genders,ones)
        return(peddat)
}

co.peddat=function(){ # cousins
        pedsize <- 8
        pedids <- rep(1,pedsize) 
        ids <- (1:pedsize) 
        dadids <- c(3,6,0,8,8,0,0,0)
        momids <- c(4,5,0,7,7,0,0,0)
        zeros <- rep(0,pedsize)
        genders <- c(2,2,1,2,2,1,2,1)
        ones <- rep(1,pedsize)
        peddat <- cbind(pedids,ids,dadids,momids,zeros,zeros,zeros,genders,ones)
        return(peddat)
}

user.peddat=function(userdat) { # user-defined
  if(is.null(userdat)) { 
    return(NULL) 
  } else {
    return(cbind(1,userdat$ids,userdat$dadids,userdat$momids,0,0,0,
                 userdat$gender,1))
  }
}

apply.errs<-function(x,hom2hom.err) {
  # x is a vector of raw's. 
  x.orig<-x
  # homozygotes called heterozygous
  x[x.orig==as.raw(1)]<-hom1.err(sum(x.orig==1),hom2hom.err)
  x[x.orig==as.raw(2)]<-het.err(sum(x.orig==2))
  x[x.orig==as.raw(3)]<-hom3.err(sum(x.orig==3),hom2hom.err)
  return(x)
}


hom1.err<-function(n,hom2hom.err) {
  # homoz 01 have prob 1-hom2hom.err of being 02 and hom2hom.err of being 03
  return(as.raw(2+rbinom(n,size=1,prob=(1-hom2hom.err))))
}

het.err<-function(n) {
  # heterozygotes equally likely to be one of two homozygotes
  return(as.raw(1+2*rbinom(n,size=1,prob=1/2)))
}

hom3.err<-function(n,hom2hom.err) {
  # homoz 03 have prob hom2hom.err of being 01 and 1-hom2hom.err of being 02
  return(as.raw(1+rbinom(n,size=1,prob=hom2hom.err)))
}

  
findLDfiles<-function(gdat) {
  return(unlist(lapply(gdat,getLDfile)))
}

getLDfile<-function(gdat) {
  return(gdat$GDparfile)
}


