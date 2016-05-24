genomic.clines <-
function(introgress.data=NULL, hi.index=NULL, loci.data=NULL,
                         sig.test=FALSE, method="permutation", n.reps=1000, classification=FALSE,
                         het.cor=TRUE, loci.touse=NULL, ind.touse=NULL){
  if (is.null(introgress.data)==TRUE | is.null(hi.index)==TRUE |
      is.null(loci.data)==TRUE) stop("error, missing input file(s)")
  ## makes sure enough reps (50) were included for generating 95% CI
  if (n.reps < 50) stop ("at least 50 reps are required for significance testing")
  ## let user know genomic.clines is working
  cat("genomic.clines is working; be patient, as this may take a while \n")
  ## this breaks up the individual count matrix
  individual.data<-NULL
  combos.touse<-NULL
  if (is.list(introgress.data)==TRUE){
    individual.data<-introgress.data[[1]]
    combos.touse<-introgress.data[[3]]
    alleles<-introgress.data[[6]]
    count.matrix<-introgress.data[[2]]
    rownames(count.matrix)<-loci.data[,1]
    if (is.null(individual.data)==FALSE){
      if (is.na(individual.data[1,2])==FALSE) colnames(count.matrix)<-individual.data[,2]
    }
  }
  else count.matrix<-introgress.data
  
  if (is.null(loci.touse)==FALSE){
    if (is.character(loci.touse)==TRUE & is.character(rownames(count.matrix))==FALSE)
      stop ("loci names were not supplied for subsetting")
  }
  else loci.touse<-1:dim(count.matrix)[1]
  
  if (is.null(ind.touse)==FALSE){
    if (is.character(ind.touse)==TRUE & is.character(colnames(count.matrix))==FALSE)
      stop ("individual names were not supplied for subsetting")
  }
  else ind.touse<-1:dim(count.matrix)[2]
  
  ## this retrieves the point estimate of hybrid index if hi.index is a data.frame
  if (is.data.frame(hi.index)==TRUE) hi.index<-hi.index[,2]
  ##subset files
  rownames(loci.data)<-rownames(count.matrix)
  temp.ind<-cbind(individual.data,hi.index)
  rownames(temp.ind)<-colnames(count.matrix)
  if (is.null(combos.touse)==FALSE){
    rownames(combos.touse)<-rownames(count.matrix)
    combos.touse<-combos.touse[loci.touse,]
  }	
  temp.ind<-temp.ind[ind.touse,]
  count.matrix<-count.matrix[loci.touse,ind.touse]
  loci.data<-loci.data[loci.touse,]
  if (is.matrix(temp.ind)==TRUE) individual.data<-temp.ind[,1:2]
  if (is.matrix(temp.ind)==TRUE) hi.index<-as.numeric(temp.ind[,3])
  else hi.index<-as.numeric(temp.ind)
  temp.ind<-NULL
  n.loci<-length(loci.touse)
  n.ind<-length(ind.touse)
  ## this code conducts multinomial regressions for the observed data (y) and hybrid index vector (x)
  if (sig.test==FALSE){
    ## data objects to save results
    AA.slope<-numeric(n.loci)
    Aa.slope<-numeric(n.loci)
    AA.int<-numeric(n.loci)
    Aa.int<-numeric(n.loci)
    AA.real.fitted.array<-array(dim=c(n.loci,n.ind))  
    Aa.real.fitted.array<-array(dim=c(n.loci,n.ind))
    aa.real.fitted.array<-array(dim=c(n.loci,n.ind))
    
    ## Information on major storage variables: 
    ## AA.real.fitted.array, Aa.real.fitted.array, and aa.real.fitted.array are two dimensional arrays (loci, individuals)
    ## that store the fitted regression values for each individual at each locus based on the observed data. Note, AA corresponds to the 2 genotype.
    ## Similarly AA(Aa, or aa).fitted.all are two dimensional arrays (reps, individuals) used to calculate the mean neutral expectations. They do not
    ## store data from multiple loci. AA(Aa, or aa).fitted.array are three dimensional arrays (loci, individuals, reps) that store fitted values for each
    ## of the neutral replicates used for significance testing. AA(Aa, or aa).neutral.ub(lb) are two dimensional arrays (loci, individuals) that store
    ## the lower and upper bounds for the 95% CI of expected genotype frequencies as a function of hybrid index for the neutral replicates. These
    ## are used for making cline plots.
    for (i in 1:n.loci){
      cat("estimating genomic cline for:",as.character(loci.data[i,1]),"\n")
      local.cnt<-count.matrix[i,]
      local.cnt<-local.cnt[!is.na(local.cnt)]
      if (length(unique(local.cnt))==1){
        warning ("warning, invariant locus included: ", as.character(loci.data[i,1]))
        clines.out<-fit.invariant.clines(count.matrix[i,],n.ind,loci.data[i,2])
        if (loci.data[i,2]=="C" | loci.data[i,2]=="c"){
          AA.real.fitted.array[i,]<-clines.out[1,]
          Aa.real.fitted.array[i,]<-clines.out[2,]
          aa.real.fitted.array[i,]<-clines.out[3,] 
        }
        else{
          AA.real.fitted.array[i,]<-clines.out[1,]
          aa.real.fitted.array[i,]<-clines.out[2,]
        }
      }
      else{
        reg.out<-multinom(count.matrix[i,]~hi.index, trace=FALSE)
        ## for dominant data
        if (loci.data[i,2]=="D" | loci.data[i,2]=="d" | loci.data[i,2]=="H" | loci.data[i,2]=="h"){
          AA.slope[i]<-coef(reg.out)[2]
          AA.int[i]<-coef(reg.out)[1]
          Hx<-exp(AA.int[i]+AA.slope[i]*hi.index)
          AA.real.fitted.array[i,]<-exp(AA.int[i]+AA.slope[i]*hi.index)/(1+Hx)
          aa.real.fitted.array[i,]<-1-AA.real.fitted.array[i,]	
        }
        ## for co-dominant data
        ## this now uses the fit.c.clines function
        else if (loci.data[i,2]=="C" | loci.data[i,2]=="c"){
          clines.out<-fit.c.clines(reg.out,hi.index,count.matrix[i,],n.ind)
          AA.real.fitted.array[i,]<-clines.out[1,]
          Aa.real.fitted.array[i,]<-clines.out[2,]
          aa.real.fitted.array[i,]<-clines.out[3,]
        }
      }
    }
    rownames(AA.real.fitted.array)<-loci.data[,1]
    rownames(Aa.real.fitted.array)<-loci.data[,1]
    rownames(aa.real.fitted.array)<-loci.data[,1]
    if (is.null(individual.data)==FALSE){
      colnames(AA.real.fitted.array)<-individual.data[,2]
      colnames(Aa.real.fitted.array)<-individual.data[,2]
      colnames(aa.real.fitted.array)<-individual.data[,2]
    }
    ## clines.out is a list, the first element gives the locus names
    ## etc., the next three elements are the fitted values for AA, Aa,
    ## and aa; the next three elements are the upper and lower bounds
    ## (2.5 and 97.5) for the simulated/neutral fitted values for AA,
    ## Aa, and aa--these are NULL if sig.test=FALSE
    summary.data<-cbind(loci.data[,1],loci.data[,2])
    colnames(summary.data)<-c("locus","type")
    clines.out<-list(summary.data, AA.real.fitted.array, Aa.real.fitted.array,
                     aa.real.fitted.array,NULL,NULL,NULL,count.matrix,hi.index,loci.data)
    names(clines.out)<-c("Summary.data","Fitted.AA","Fitted.Aa","Fitted.aa","Neutral.AA","Neutral.Aa",
                         "Neutral.aa","Count.matrix","hybrid.index","Loci.data")
    cat ("genomic clines analysis complete! \n")
    return(clines.out)
  }
  ## significance test for cline analysis
  else if (sig.test==TRUE){
    ## verify that a proper method was specified
    if (method!="permutation" & method!="parametric") stop("invalid method specified")
    cat ("begining analysis. \n")
    ## set up data objects for results
    p.value<-numeric(n.loci)
    lnlik.ratios<-numeric(n.loci)
    AA.slope<-numeric(n.loci)
    Aa.slope<-numeric(n.loci)
    AA.int<-numeric(n.loci)
    Aa.int<-numeric(n.loci)
    
    AA.fitted.array<-array(rep(NA,n.ind*n.reps),dim=c(n.ind,n.reps))
    Aa.fitted.array<-array(rep(NA,n.ind*n.reps),dim=c(n.ind,n.reps))
    aa.fitted.array<-array(rep(NA,n.ind*n.reps),dim=c(n.ind,n.reps))
    
    AA.real.fitted.array<-array(rep(NA,n.loci*n.ind),dim=c(n.loci,n.ind))
    Aa.real.fitted.array<-array(rep(NA,n.loci*n.ind),dim=c(n.loci,n.ind))
    aa.real.fitted.array<-array(rep(NA,n.loci*n.ind),dim=c(n.loci,n.ind))
    
    AA.neutral.ub<-array(rep(NA,n.loci*n.ind),dim=c(n.loci,n.ind))
    Aa.neutral.ub<-array(rep(NA,n.loci*n.ind),dim=c(n.loci,n.ind))
    aa.neutral.ub<-array(rep(NA,n.loci*n.ind),dim=c(n.loci,n.ind))

    AA.neutral.lb<-array(rep(NA,n.loci*n.ind),dim=c(n.loci,n.ind))
    Aa.neutral.lb<-array(rep(NA,n.loci*n.ind),dim=c(n.loci,n.ind))
    aa.neutral.lb<-array(rep(NA,n.loci*n.ind),dim=c(n.loci,n.ind))
    ## save genotype probabilities
    AA.realProbQ<-numeric(n.loci)
    Aa.realProbQ<-numeric(n.loci)
    aa.realProbQ<-numeric(n.loci)
    
    ## permutation neutral simulations
    if (method=="permutation"){
      ## simulation of large population under neutrality
      sam.neutral<-array(dim=c(n.ind,n.reps))
      for (j in 1:n.ind){
        for(k in 1:n.reps){
          sam.neutral[j,k]<-sample(count.matrix[,j],1,replace=FALSE)
        }	
      }
      ## simulation of neutral replicate populations
      sam.gen<-array(dim=c(n.ind,n.reps))
      for (j in 1:n.ind){
        for(k in 1:n.reps){
          sam.gen[j,k]<-sample(count.matrix[,j],1,replace=FALSE)
        }	
      }			
    }	
    ## parametric neutral simulations
    if (method=="parametric"){
      ## set up data objects for results
      sam.neutral<-array(dim=c(n.ind,n.reps))
      sam.gen<-array(dim=c(n.ind,n.reps))
      
      ## E.A is the prob of sampling an allele from allelic class 2,
      ## which is at higher frequency in spa (parental1, species with
      ## low hybrid index)
      
      E.A<-array(rep(NA,n.loci*n.ind),dim=c(n.loci,n.ind))
      if (is.null(combos.touse)==FALSE){
        for (i in 1:n.loci){
          ## prior to 16 Oct 09
          ## E.A[i,]<-combos.touse$spb.c1[i] + (combos.touse$spa.c1[i] - combos.touse$spb.c1[i]) * hi.index
          ## now
          E.A[i,]<-combos.touse$spa.c2[i] + (combos.touse$spb.c2[i] - combos.touse$spa.c2[i]) * hi.index
        }
      }			
      else if (is.null(combos.touse)==TRUE){
        ## tests for combos.touse, if it is still NULL fixed
        ## differences are assumed, with spb (high hybrid index)
        ## having a frequency of 0 and spa (low hybrid index) having a
        ## frequency of 1
        for (i in 1:n.loci){
          E.A[i,]<-1-hi.index
        }
      }
      ## calculate excess het if het.cor==TRUE and at least 10 co-dominant markers are included
      if (het.cor==TRUE & (sum(loci.data[,2]=="C") + (sum(loci.data[,2]=="c"))) >=10){
        co.markers<-which(loci.data[,2]=="C" | loci.data[,2]=="c")
        obs.het<-numeric(n.ind)
        for(j in 1:n.ind){
          obs.het[j]<-sum(count.matrix[co.markers,j]==1,na.rm=TRUE)/
            sum(is.na(count.matrix[co.markers,j])==FALSE)
        }
        exp.het<-numeric(n.ind)
        for (j in 1:n.ind){
          exp.het[j]<-mean(2*E.A[co.markers,j]*(1-E.A[co.markers,j]))
        }
        excess.het<-obs.het-exp.het
      }
      else  excess.het<-0
      ## set up E.AA, E.Aa, and E.aa arrays
      E.AA<-array(dim=c(n.loci,n.ind))
      E.aa<-array(dim=c(n.loci,n.ind))
      if ((sum(loci.data[,2]=="C") + (sum(loci.data[,2]=="c"))) > 0)
        E.Aa<-array(dim=c(n.loci,n.ind))
      
      ## calculate E.AA, E.Aa, E.aa for co-dominant markers, or for
      ## dominant and haploid markers E.AA, E.a
      for (i in 1:n.loci){
        if (loci.data[i,2]=="D" | loci.data[i,2]=="d"){
          ## determines if presence allele is present at higher
          ## freq. in parental1 (spa), if so E.A is presence allele
          if ((alleles[i,1]=="1" & combos.touse[i,3]=="2") |
              (alleles[i,1]=="0" & combos.touse[i,3]=="1")){
            E.AA[i,]<-E.A[i,]^2 + 2*E.A[i,]*(1-E.A[i,])
            E.aa[i,]<-(1-E.A[i,])^2
          }
          ## determines if absence allele is present at higher
          ## freq. in parental1 (spa)
          else if ((alleles[i,1]=="1" & combos.touse[i,3]=="1") |
                   (alleles[i,1]=="0" & combos.touse[i,3]=="2")){
            E.AA[i,]<-E.A[i,]^2
            E.aa[i,]<-(1-E.A[i,])^2 + 2*E.A[i,]*(1-E.A[i,])
          }
          else
            stop(paste("dominant marker data is not coded properly for marker:",
                       loci.data[i,1]))
        }	
        else if (loci.data[i,2]=="H" | loci.data[i,2]=="h"){
          E.AA[i,]<-E.A[i,]
          E.aa[i,]<-1-E.A[i,]
        }	
        else if	(loci.data[i,2]=="C" | loci.data[i,2]=="c"){
          ## calculate E.AA etc. and normalize
          E.AA[i,]<-E.A[i,]^2 - (E.A[i,]^2/(E.A[i,]^2 + (1-E.A[i,])^2))*excess.het
          E.Aa[i,]<-2*E.A[i,]*(1-E.A[i,]) + excess.het
          E.aa[i,]<-(1-E.A[i,])^2 - ((1-E.A[i,])^2/(E.A[i,]^2 + (1-E.A[i,])^2))*excess.het
          for(j in 1:n.ind){
            if(E.aa[i,j]<0) E.aa[i,j]<-0
            else if(E.aa[i,j]>1) E.aa[i,j]<-1

            if(E.AA[i,j]<0) E.AA[i,j]<-0
            else if(E.AA[i,j]>1) E.AA[i,j]<-1

            if(E.Aa[i,j]<0) E.Aa[i,j]<-0
            else if(E.Aa[i,j]>1) E.Aa[i,j]<-1
          }	
        }	
      }
    }	

    ##begin cycling through loci			
    for (i in 1:n.loci){
      cat("estimating genomic cline for:",as.character(loci.data[i,1]),"\n")
      ## locus specific neutral samples for parametric method
      if (method=="parametric"){
        ## clear sams
        sam.neutral[,]<-NA
        sam.gen[,]<-NA
        if (loci.data[i,2]=="D" | loci.data[i,2]=="d"){
          ## simulation of large population under neutrality
          for (j in 1:n.ind){
            for(k in 1:n.reps){
              sam.neutral[j,k]<-sample(c(1,0),1,replace=TRUE,prob=c(E.AA[i,j],E.aa[i,j]))
            }	
          }
          ## simulation of neutral replicate populations
          for (j in 1:n.ind){
            for(k in 1:n.reps){
              sam.gen[j,k]<-sample(c(1,0),1,replace=TRUE,prob=c(E.AA[i,j],E.aa[i,j]))
            }	
          }         
        }
        else if (loci.data[i,2]=="H" | loci.data[i,2]=="h"){
          ## simulation of large population under neutrality
          for (j in 1:n.ind){
            for(k in 1:n.reps){
              sam.neutral[j,k]<-sample(c(1,0),1,replace=TRUE,prob=c(E.AA[i,j],E.aa[i,j]))
            }	
          }
          ## simulation of neutral replicate populations
          for (j in 1:n.ind){
            for(k in 1:n.reps){
              sam.gen[j,k]<-sample(c(1,0),1,replace=TRUE,prob=c(E.AA[i,j],E.aa[i,j]))
            }	
          }
        }
        else if	(loci.data[i,2]=="C" | loci.data[i,2]=="c"){
          ## simulation of large population under neutrality
          for (j in 1:n.ind){
            for(k in 1:n.reps){
              sam.neutral[j,k]<-sample(c(2,1,0),1,replace=TRUE,
                                         prob=c(E.AA[i,j],E.Aa[i,j],E.aa[i,j]))
            }	
          }
          ##simulation of neutral replicate populations
          for (j in 1:n.ind){
            for(k in 1:n.reps){
              sam.gen[j,k]<-sample(c(2,1,0),1,replace=TRUE,
                                     prob=c(E.AA[i,j],E.Aa[i,j],E.aa[i,j]))
            }	
          }
        }
      }
      ## clear fitted.array
      AA.fitted.array[,]<-NA
      Aa.fitted.array[,]<-NA
      aa.fitted.array[,]<-NA
      ## clines analysis for mean neutral expectations
      AA.fitted.all<-array(dim=c(n.reps,n.ind))
      Aa.fitted.all<-array(dim=c(n.reps,n.ind))
      aa.fitted.all<-array(dim=c(n.reps,n.ind))
      for (k in 1:n.reps){
        local.cnt<-sam.neutral[,k]
        local.cnt<-local.cnt[!is.na(local.cnt)]
        if(length(unique(local.cnt))==1){
          clines.out<-fit.invariant.clines(sam.neutral[,k],n.ind,loci.data[i,2])
          if (loci.data[i,2]=="C" | loci.data[i,2]=="c"){
            AA.fitted.all[k,]<-clines.out[1,]
            Aa.fitted.all[k,]<-clines.out[2,]
            aa.fitted.all[k,]<-clines.out[3,]

          }
          else{
            AA.fitted.all[k,]<-clines.out[1,]
            aa.fitted.all[k,]<-clines.out[2,]
          }
        }
        else {
          reg.out<-multinom(sam.neutral[,k]~hi.index, trace=FALSE)
          ## for dominant or haploid markers data
          if (loci.data[i,2]=="D" | loci.data[i,2]=="d" | loci.data[i,2]=="H" |
              loci.data[i,2]=="h"){
            AA.slope.neutral<-coef(reg.out)[2]
            AA.int.neutral<-coef(reg.out)[1]
            Hx<-exp(AA.int.neutral+AA.slope.neutral*hi.index)
            AA.fitted.all[k,]<-exp(AA.int.neutral+AA.slope.neutral*hi.index)/(1+Hx)
            aa.fitted.all[k,]<-1-AA.fitted.all[k,]
          }	
          ## for co-dominant data
          ## this now uses the fit.c.clines function
          else if (loci.data[i,2]=="C" | loci.data[i,2]=="c"){
            clines.out<-fit.c.clines(reg.out,hi.index,sam.neutral[,k],n.ind)
            AA.fitted.all[k,]<-clines.out[1,]
            Aa.fitted.all[k,]<-clines.out[2,]
            aa.fitted.all[k,]<-clines.out[3,]
          }
        }
      }
      ## sum results from each rep
      if (loci.data[i,2]=="D" | loci.data[i,2]=="d" | loci.data [i,2]=="H" |
          loci.data [i,2]=="h"){
        AA.fitted.mean<-apply(AA.fitted.all,2,mean,na.rm=TRUE)
        aa.fitted.mean<-apply(aa.fitted.all,2,mean,na.rm=TRUE)
      }	
      else if (loci.data[i,2]=="C" | loci.data[i,2]=="c"){
        AA.fitted.mean<-apply(AA.fitted.all,2,mean,na.rm=TRUE)
        Aa.fitted.mean<-apply(Aa.fitted.all,2,mean,na.rm=TRUE)
        aa.fitted.mean<-apply(aa.fitted.all,2,mean,na.rm=TRUE)
      }
      ## reset slope/int
      AA.slope.neutral<-NA
      AA.int.neutral<-NA
      Aa.slope.neutral<-NA
      Aa.int.neutral<-NA
      
      ## cline analysis for individual neutral expectations incorporating sampling error

      ## set up likelihood data objects
      prob.obs.model1<-rep(NA,n.ind)
      prob.obs.model0<-rep(NA,n.ind)
      ln.likelihood.model1<-numeric(n.reps)
      ln.likelihood.model0<-numeric(n.reps)
      ln.likelihood.ratio10<-numeric(n.reps)

      ## set up data objects to save genotype specific probabilities
      AA.prob<-numeric(n.reps)
      Aa.prob<-numeric(n.reps)
      aa.prob<-numeric(n.reps)
      
      ## perform regressions			
      for (k in 1:n.reps){
        local.cnt<-sam.gen[,k]
        local.cnt<-local.cnt[!is.na(local.cnt)]
        if(length(unique(local.cnt))==1){
          clines.out<-fit.invariant.clines(sam.gen[,k],n.ind,loci.data[i,2])
          if (loci.data[i,2]=="C" | loci.data[i,2]=="c"){
            AA.fitted.array[,k]<-clines.out[1,]
            Aa.fitted.array[,k]<-clines.out[2,]
            aa.fitted.array[,k]<-clines.out[3,]
            ## save genotype probabilities
            if(classification==TRUE){
              AA.prob[k]<-sum(AA.fitted.array[,k],na.rm=TRUE)
              Aa.prob[k]<-sum(Aa.fitted.array[,k],na.rm=TRUE)
              aa.prob[k]<-sum(aa.fitted.array[,k],na.rm=TRUE)
            }
          }
          else{
            AA.fitted.array[,k]<-clines.out[1,]
            aa.fitted.array[,k]<-clines.out[2,]
          }
        }
        else{
          reg.out<-multinom(sam.gen[,k]~hi.index, trace=FALSE)
          ## for dominant or haploid marker data
          if (loci.data[i,2]=="D" | loci.data[i,2]=="d" | loci.data[i,2]=="H" | loci.data[i,2]=="h"){
            AA.slope.neutral<-coef(reg.out)[2]
            AA.int.neutral<-coef(reg.out)[1]
            Hx<-exp(AA.int.neutral+AA.slope.neutral*hi.index)
            AA.fitted.array[,k]<-exp(AA.int.neutral+AA.slope.neutral*hi.index)/(1+Hx)
            aa.fitted.array[,k]<-1-AA.fitted.array[,k]
            ## save genotype probabilities
            if(classification==TRUE){
              AA.prob[k]<-sum(AA.fitted.array[,k],na.rm=TRUE)
              aa.prob[k]<-sum(aa.fitted.array[,k],na.rm=TRUE)
            }
          }
          ## for co-dominant data
          ## this now uses the fit.c.clines function
          else if (loci.data[i,2]=="C" | loci.data[i,2]=="c"){
            clines.out<-fit.c.clines(reg.out,hi.index,sam.gen[,k],n.ind)
            AA.fitted.array[,k]<-clines.out[1,]
            Aa.fitted.array[,k]<-clines.out[2,]
            aa.fitted.array[,k]<-clines.out[3,]
            ## save genotype probabilities
            if(classification==TRUE){
              AA.prob[k]<-sum(AA.fitted.array[,k],na.rm=TRUE)
              Aa.prob[k]<-sum(Aa.fitted.array[,k],na.rm=TRUE)
              aa.prob[k]<-sum(aa.fitted.array[,k],na.rm=TRUE)
            }
          }
        }
        ## calculates probability of models
        ## for dominant or haploid marker data
        if (loci.data[i,2]=="D" | loci.data[i,2]=="d" | loci.data[i,2]=="H" |
            loci.data[i,2]=="h"){
          for(z in 1:n.ind){
            if(is.na(count.matrix[i,z])==FALSE){
              if(is.na(sam.gen[z,k])==FALSE){
                if(sam.gen[z,k]==2) {
                  prob.obs.model1[z]<-AA.fitted.array[z,k]
                  prob.obs.model0[z]<-AA.fitted.mean[z]
                }	
                else if(sam.gen[z,k]==0){ 
                  prob.obs.model1[z]<-aa.fitted.array[z,k]
                  prob.obs.model0[z]<-aa.fitted.mean[z]
                }	
              }
              else if(is.na(sam.gen[z,k])==TRUE){
                prob.obs.model1[z]<-0.50
                prob.obs.model0[z]<-0.50
              }	
            }	
          }
        }
        ## for co-dominant data
        if (loci.data[i,2]=="C" | loci.data[i,2]=="c"){
          for(z in 1:n.ind){
            if(is.na(count.matrix[i,z])==FALSE){
              if(is.na(sam.gen[z,k])==FALSE){
                if(sam.gen[z,k]==2) {
                  prob.obs.model1[z]<-AA.fitted.array[z,k]
                  prob.obs.model0[z]<-AA.fitted.mean[z]
                }	
                else if(sam.gen[z,k]==1) {
                  prob.obs.model1[z]<-Aa.fitted.array[z,k]
                  prob.obs.model0[z]<-Aa.fitted.mean[z]
                }	
                else if(sam.gen[z,k]==0){ 
                  prob.obs.model1[z]<-aa.fitted.array[z,k]
                  prob.obs.model0[z]<-aa.fitted.mean[z]
                }	
              }
              else if(is.na(sam.gen[z,k])==TRUE){
                prob.obs.model1[z]<-0.33
                prob.obs.model0[z]<-0.33
              }	
            }	
          }
        }		
        ln.likelihood.model1[k]<-sum(log(prob.obs.model1),na.rm=TRUE)
        ln.likelihood.model0[k]<-sum(log(prob.obs.model0),na.rm=TRUE)
        ln.likelihood.ratio10[k]<-ln.likelihood.model1[k]-ln.likelihood.model0[k]
      }
      ## save upper and lower bounds of 95%CI from neutral
      ## simulations, these are for making cline plots
      for (j in 1:n.ind){
        AA.neutral.ub[i,j]<-sort(AA.fitted.array[j,])[n.reps*0.975]
        Aa.neutral.ub[i,j]<-sort(Aa.fitted.array[j,])[n.reps*0.975]
        aa.neutral.ub[i,j]<-sort(aa.fitted.array[j,])[n.reps*0.975]
        AA.neutral.lb[i,j]<-sort(AA.fitted.array[j,])[n.reps*0.0275]
        Aa.neutral.lb[i,j]<-sort(Aa.fitted.array[j,])[n.reps*0.0275]
        aa.neutral.lb[i,j]<-sort(aa.fitted.array[j,])[n.reps*0.0275]
      }
      ## assign rownames, this is necessary for indexing by locus names 
      ## when making cline plots
      rownames(AA.neutral.ub)<-loci.data[,1]
      rownames(Aa.neutral.ub)<-loci.data[,1]
      rownames(aa.neutral.ub)<-loci.data[,1]
      rownames(AA.neutral.lb)<-loci.data[,1]
      rownames(Aa.neutral.lb)<-loci.data[,1]
      rownames(aa.neutral.lb)<-loci.data[,1]
      local.cnt<-count.matrix[i,]
      local.cnt<-local.cnt[!is.na(local.cnt)]
      if(length(unique(local.cnt))==1){
        warning ("warning, invariant locus included: ", as.character(loci.data[i,1]))
        clines.out<-fit.invariant.clines(count.matrix[i,],n.ind,loci.data[i,2])
        if (loci.data[i,2]=="C" | loci.data[i,2]=="c"){
          AA.real.fitted.array[i,]<-clines.out[1,]
          Aa.real.fitted.array[i,]<-clines.out[2,]
          aa.real.fitted.array[i,]<-clines.out[3,]
          if(classification==TRUE){
            AA.realProb<-sum(AA.real.fitted.array[i,],na.rm=TRUE)
            Aa.realProb<-sum(Aa.real.fitted.array[i,],na.rm=TRUE)
            aa.realProb<-sum(aa.real.fitted.array[i,],na.rm=TRUE)
          }
        }
        else{
          AA.real.fitted.array[i,]<-clines.out[1,]
          aa.real.fitted.array[i,]<-clines.out[2,]
          if(classification==TRUE){
            AA.realProb<-sum(AA.real.fitted.array[i,],na.rm=TRUE)
            aa.realProb<-sum(aa.real.fitted.array[i,],na.rm=TRUE) 
          }
        }
      }
      else{
        ## multinomial regression on observed data
        reg.out<-multinom(count.matrix[i,]~hi.index, trace=FALSE)
        ## for dominant or haploid marker data
        if (loci.data[i,2]=="D" | loci.data[i,2]=="d" | loci.data[i,2]=="H" |
            loci.data[i,2]=="h"){
          AA.slope[i]<-coef(reg.out)[2]
          AA.int[i]<-coef(reg.out)[1]
          Hx<-exp(AA.int[i]+AA.slope[i]*hi.index)
          AA.real.fitted.array[i,]<-exp(AA.int[i]+AA.slope[i]*hi.index)/(1+Hx)
          aa.real.fitted.array[i,]<-1-AA.real.fitted.array[i,]
          ## save genotype probabilities
          if(classification==TRUE){
            AA.realProb<-sum(AA.real.fitted.array[i,],na.rm=TRUE)
            aa.realProb<-sum(aa.real.fitted.array[i,],na.rm=TRUE)
          }
        }	
        ## for co-dominant data
        else if (loci.data[i,2]=="C" | loci.data[i,2]=="c"){
          clines.out<-fit.c.clines(reg.out,hi.index,count.matrix[i,],n.ind)
          AA.real.fitted.array[i,]<-clines.out[1,]
          Aa.real.fitted.array[i,]<-clines.out[2,]
          aa.real.fitted.array[i,]<-clines.out[3,]
          ## save genotype probabilities
          if(classification==TRUE){
            AA.realProb<-sum(AA.real.fitted.array[i,],na.rm=TRUE)
            Aa.realProb<-sum(Aa.real.fitted.array[i,],na.rm=TRUE)
            aa.realProb<-sum(aa.real.fitted.array[i,],na.rm=TRUE)
          }
        }
      }
      ## calculate probabilities
      prob.obs.model1<-(rep(NA,n.ind))
      prob.obs.model0<-(rep(NA,n.ind))
      if (loci.data[i,2]=="D" | loci.data[i,2]=="d" | loci.data[i,2]=="H" |
          loci.data[i,2]=="h"){
        for(z in 1:n.ind){
          if(is.na(count.matrix[i,z])==FALSE){
            if(count.matrix[i,z]==1) {
              prob.obs.model1[z]<-AA.real.fitted.array[i,z]
              prob.obs.model0[z]<-AA.fitted.mean[z]
            }	
            else if(count.matrix[i,z]==0) {
              prob.obs.model1[z]<-aa.real.fitted.array[i,z]
              prob.obs.model0[z]<-aa.fitted.mean[z]
            }	
          }
        }
      }			
      if (loci.data[i,2]=="C" | loci.data[i,2]=="c"){
        for(z in 1:n.ind){
          if(is.na(count.matrix[i,z])==FALSE){
            if(count.matrix[i,z]==2) {
              prob.obs.model1[z]<-AA.real.fitted.array[i,z]
              prob.obs.model0[z]<-AA.fitted.mean[z]
            }	
            else if(count.matrix[i,z]==1) {
              prob.obs.model1[z]<-Aa.real.fitted.array[i,z]
              prob.obs.model0[z]<-Aa.fitted.mean[z]
            }	
            else if(count.matrix[i,z]==0) {
              prob.obs.model1[z]<-aa.real.fitted.array[i,z]
              prob.obs.model0[z]<-aa.fitted.mean[z]
            }	
          }	
        }
      }	
      ln.likelihood.model1.real<-sum(log(prob.obs.model1),na.rm=TRUE)
      ln.likelihood.model0.real<-sum(log(prob.obs.model0),na.rm=TRUE)
      ln.likelihood.ratio10.real<-ln.likelihood.model1.real-ln.likelihood.model0.real
      ## calculate and save p.value and ln.lik ratio	
      p.value[i]<-(sum(ln.likelihood.ratio10>=ln.likelihood.ratio10.real))/n.reps
      lnlik.ratios[i]<-ln.likelihood.ratio10.real
      ## calculate and save genotype specific quantiles
      if(classification==TRUE){
        AA.realProbQ[i]<-sum(AA.prob < AA.realProb, na.rm=TRUE)/n.reps
        if (loci.data[i,2]=="C" | loci.data[i,2]=="c"){
          Aa.realProbQ[i]<-sum(Aa.prob < Aa.realProb, na.rm=TRUE)/n.reps
          aa.realProbQ[i]<-sum(aa.prob < aa.realProb, na.rm=TRUE)/n.reps
        }
      }
      ##the bracket below closes the locus loop
    }
  }
  ## make list for output: clines.out is a list, the first element
  ## gives the locus names, the next three elements are the fitted
  ## values for AA, Aa, and aa; the next three elements are the
  ## upper and lower bounds (2.5 and 97.5) for the simulated/neutral
  ## fitted values for AA, Aa, and aa--these are NULL if
  ## sig.test=FALSE
  AA.bounds<-list(AA.neutral.ub,AA.neutral.lb)
  Aa.bounds<-list(Aa.neutral.ub,Aa.neutral.lb)
  aa.bounds<-list(aa.neutral.ub,aa.neutral.lb)
  
  ## name columns and rows and build list object
  rownames(AA.real.fitted.array)<-loci.data[,1]
  rownames(Aa.real.fitted.array)<-loci.data[,1]
  rownames(aa.real.fitted.array)<-loci.data[,1]
  if (is.null(individual.data)==FALSE){
    colnames(AA.real.fitted.array)<-individual.data[,2]
    colnames(Aa.real.fitted.array)<-individual.data[,2]
    colnames(aa.real.fitted.array)<-individual.data[,2]
  }
  summary.data<-cbind(loci.data[,1],loci.data[,2],lnlik.ratios,p.value)
  if(classification==TRUE){
    quantiles<-cbind(AA.realProbQ,Aa.realProbQ,aa.realProbQ)
    rownames(quantiles)<-loci.data[,1]
  }
  else quantiles<-NULL
  colnames(summary.data)<-c("locus","type","lnL ratio","P-value")
  clines.out<-list(summary.data, AA.real.fitted.array, Aa.real.fitted.array,
                   aa.real.fitted.array,AA.bounds,Aa.bounds,aa.bounds,count.matrix,
                   hi.index,loci.data,quantiles)
  names(clines.out)<-c("Summary.data","Fitted.AA","Fitted.Aa","Fitted.aa",
                       "Neutral.AA","Neutral.Aa","Neutral.aa",
                       "Count.matrix","hybrid.index","Loci.data","Quantiles")
  cat ("genomic clines analysis complete! \n")
  return(clines.out)
}

