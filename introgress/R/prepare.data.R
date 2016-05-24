prepare.data <-
function(admix.gen=NULL,loci.data=NULL,
                       parental1=NULL, parental2=NULL, pop.id=TRUE,
                       ind.id=TRUE,fixed=FALSE, sep.rows=FALSE,
                       sep.columns=FALSE){
  ##the genetics library is needed for this function
  ##library(genetics)
  ## This makes sure that the genotype matrix, locus data, and parental
  ## population data were supplied
  if (is.null(dim(admix.gen))==TRUE)
    stop("Genotype matrix was not supplied")
  else if (is.null(dim(loci.data))==TRUE)
    stop("Locus information was not supplied")
  else if (is.null(parental1)==TRUE | is.null(parental2)==TRUE)
    stop("Parental population information was not supplied")
  if (is.null(colnames(loci.data))) stop("column names for loci.data are missing.")
  ## Let user know analysis is working
  cat("prepare.data is working; this may take a moment", fill=TRUE)
  ## population and individual data are stored as a matrix and removed from admix.gen
  admix.gen<-as.matrix(admix.gen)
  if (sep.columns==TRUE & sep.rows==TRUE)
    stop("sep.columns and sep.rows cannot both be true")

  if (pop.id==TRUE & ind.id==TRUE){
    if(sep.columns==TRUE){
      pop<-as.character(admix.gen[1, seq(1, ncol(admix.gen), 2)])
      ind<-as.character(admix.gen[2, seq(1, ncol(admix.gen), 2)])
    }
    else{
      pop<-as.character(admix.gen[1,])
      ind<-as.character(admix.gen[2,])
    }
    admix.gen<-as.matrix(admix.gen[3:nrow(admix.gen),])
    individual.data<-cbind(pop,ind)
  }	
  else if (pop.id==TRUE & ind.id==FALSE){
    if(sep.columns==TRUE){
      pop<-as.character(admix.gen[1, seq(1, ncol(admix.gen), 2)])
    }
    else{
      pop<-as.character(admix.gen[1,])
    }
    admix.gen<-as.matrix(admix.gen[2:dim(admix.gen)[1],])
    individual.data<-cbind(pop,rep(NA,length(pop)))
  }	
  else if (pop.id==FALSE & ind.id==TRUE){
    if(sep.columns==TRUE){
      ind<-as.character(admix.gen[1, seq(1, ncol(admix.gen), 2)])
    }
    else{
      ind<-as.character(admix.gen[1,])
    }
    admix.gen<-as.matrix(admix.gen[2:dim(admix.gen)[1],])
    individual.data<-cbind(rep(NA,length(ind)),ind)
  }	
  else individual.data<-NULL
  ## this reformats admix.gen and parentals if individual data is in
  ## two rows or columns  
  if (sep.rows==TRUE){
    ## individuals are in columns with 2 rows per locus
    admix2<-matrix(nrow=nrow(admix.gen)/2, ncol=ncol(admix.gen))
    for (i in seq(1, nrow(admix.gen), 2)){
      admix2[(i+1)/2,]<-paste(admix.gen[i,],admix.gen[i+1,],sep="/")
    }
    admix.gen<-admix2
    p1<-matrix(nrow=nrow(parental1)/2, ncol=ncol(parental1))
    for (i in seq(1, nrow(parental1), 2)){
      p1[(i+1)/2,]<-paste(parental1[i,], parental1[i+1,],sep="/")
    }
    parental1<-p1
    p2<-matrix(nrow=nrow(parental2)/2, ncol=ncol(parental2))
    for (i in seq(1, nrow(parental2), 2)){
      p2[(i+1)/2,]<-paste(parental2[i,], parental2[i+1,],sep="/")
    }
    parental2<-p2
  }
  else if (sep.columns==TRUE){
    ## individuals are in columns with 2 columns for each ind
    admix2<-matrix(nrow=nrow(admix.gen), ncol=ncol(admix.gen)/2)
    for (i in seq(1, ncol(admix.gen), 2)){
      admix2[,(i+1)/2]<-paste(admix.gen[,i],admix.gen[,i+1],sep="/")
    }
    admix.gen<-admix2
    p1<-matrix(nrow=nrow(parental1), ncol=ncol(parental1)/2)
    for (i in seq(1, ncol(parental1), 2)){
      p1[,(i+1)/2]<-paste(parental1[,i], parental1[,i+1],sep="/")
    }
    parental1<-p1
    p2<-matrix(nrow=nrow(parental2), ncol=ncol(parental2)/2)
    for (i in seq(1, ncol(parental2), 2)){
      p2[,(i+1)/2]<-paste(parental2[,i], parental2[,i+1],sep="/")
    }
    parental2<-p2
  }
  
  ## this verifies that the same number of loci are present in all data objects
  if (is.matrix(parental1)==TRUE)
    test.data.objects(admix.gen,loci.data,parental1,parental2)

  ## test for alleles not encountered in parental populations and
  ## replace with NAs
  if (is.matrix(parental1)==TRUE)
    admix.gen<-test.genotypes(admix.gen,loci.data,parental1,parental2)
  
  ## fix up file formats
  if (fixed==FALSE) parental1<-as.matrix(parental1)
  if (fixed==FALSE) parental2<-as.matrix(parental2)		

  ## determine number of loci and number of individuals
  n.loci<-dim(admix.gen)[1]
  n.ind<-dim(admix.gen)[2]

  ## make sure NA's are correct in admixed and parentals
  for(i in 1:n.loci){
    for(j in 1:n.ind){
      if(is.na(admix.gen[i,j])==FALSE){
        if(admix.gen[i,j]=="NA"| admix.gen[i,j]=="NA/NA") admix.gen[i,j]<-NA
      }
    }
  }

  if (fixed==FALSE){
    for(i in 1:n.loci){
      for(j in 1:dim(parental1)[2]){
        if(is.na(parental1[i,j])==FALSE){
          if(parental1[i,j]=="NA" | parental1[i,j]=="NA/NA") parental1[i,j]<-NA
        }
      }
    }
    for(i in 1:n.loci){
      for(j in 1:dim(parental2)[2]){
        if(is.na(parental2[i,j])==FALSE){
          if(parental2[i,j]=="NA"| parental2[i,j]=="NA/NA") parental2[i,j]<-NA
        }
      }
    }    
  }

  cat("Processing data for", n.ind, "individuals and", n.loci, "loci.", fill=TRUE)

  ## if the parental allele frequencies exhibit fixed differences
  ## "fixed" should be set to TRUE and parental1 and parental2 should
  ## each be a single value giving the  character that was used to code the
  ## allele derived from each population.  If fixed=TRUE the following
  ## code is used for computing the count.matrix
  if (fixed==TRUE){
    ## parental data must be supplied for dominant markers so that allele frequencies can be estimated
    if (sum(loci.data[,2]=="D") | sum(loci.data[,2]=="d"))
      stop("parental data must be supplied for dominant markers")

    ## make count matrix, parental1 is given a value of 2
    count.matrix<-array(dim=c(n.loci,n.ind))
    for (i in 1:n.loci){
      count.matrix[i,]<-allele.count(as.genotype(admix.gen[i,],allow.partial.missing=TRUE),
                                     parental1,na.rm=TRUE)[1:length(admix.gen[i,])]
      count.matrix[i, is.na(admix.gen[i,])]<-NA
    }
    rownames(count.matrix)<-as.character(loci.data[,1])
    if (is.null(individual.data)==FALSE){
      if (is.na(individual.data[1,2])==FALSE) colnames(count.matrix)<-as.character(individual.data[,2])
    }
    introgress.data<-list(Individual.data=individual.data, Count.matrix=count.matrix,
                          Combos.to.use=NULL, Parental1.allele.freq=NULL,
                          Parental2.allele.freq=NULL, Alleles=NULL, Admix.gen=admix.gen)
    return(introgress.data)
  }
  
  ## if the parental allele frequencies do not exhibit fixed
  ## differences fixed should be set to FALSE (the default value) and
  ## parental1 and parental2 should be a matrix with rows for loci,
  ## columns for individuals and individual genotypes for each (similar
  ## to admix.gen)

  ## if fixed=FALSE the following code is used for computing the count.matrix
  else{
    ## first parental allele frequencies are calculated for each locus for parental1 and parental2
    n.ind.p1<-dim(parental1)[2]
    p1.allele.freq<-NULL
    n.ind.p2<-dim(parental1)[2]
    p2.allele.freq<-NULL
    ploidy<-numeric(n.loci)
    for(i in 1:n.loci){
      if (loci.data[i,2]=="C" | loci.data[i,2]=="c") ploidy[i]<-2
      else if (loci.data[i,2]=="D" | loci.data[i,2]=="d" | loci.data[i,2]=="H" |
               loci.data[i,2]=="h") ploidy[i]<-1
      temp.count<-allele.count(as.genotype(c(parental1[i,],parental2[i,]),
                                           allow.partial.missing=TRUE),na.rm=TRUE)
      p1.freq<-sum(temp.count[(1:n.ind.p1),1])
      p2.freq<-sum(temp.count[(n.ind.p1+1):dim(temp.count)[1],1])
      if (dim(temp.count)[2] < 2) seq<-1
      else seq<-2:dim(temp.count)[2]
      for (j in seq){
        p1.freq2<-sum(temp.count[(1:n.ind.p1),j])
        p2.freq2<-sum(temp.count[(n.ind.p1+1):dim(temp.count)[1],j])
        p1.freq<-c(p1.freq,p1.freq2)
        p2.freq<-c(p2.freq,p2.freq2)
      }
      p1.total<-sum(p1.freq)
      p2.total<-sum(p2.freq)
      p1.freq<-p1.freq/p1.total
      p2.freq<-p2.freq/p2.total
      ## Convert from band frequencies to allele frequencies assuming HWE
      if (loci.data[i,2]=="D" | loci.data[i,2]=="d"){
        alleles<-allele.names(as.genotype(c(parental1[i,],parental2[i,]),allow.partial.missing=TRUE))
        if (alleles[1]=="0"){
          p1.freq[1]<-sqrt(p1.freq[1])
          p2.freq[1]<-sqrt(p2.freq[1])
          p1.freq[2]<-1-p1.freq[1]
          p2.freq[2]<-1-p2.freq[1]
        }
        else if (alleles[1]=="1"){
          p1.freq[2]<-sqrt(p1.freq[2])
          p2.freq[2]<-sqrt(p2.freq[2])
          p1.freq[1]<-1-p1.freq[2]
          p2.freq[1]<-1-p2.freq[2]
        }  
      }  
      ## first time through start building matrix of frequencies
      if (i==1) {
        p1.freq.all<-p1.freq
        p2.freq.all<-p2.freq
      }
      ## after first time new locus is added on as a new row, if it
      ## has more alleles than other loci, NAs must be added
      else { ## i.e. (i > 1)
        if (i==2) size<-length(p1.freq.all)
        else size<-dim(p1.freq.all)[2]
        
        ## if new locus has more alleles add NAs
        if (length(p1.freq) > size){
          dif<-length(p1.freq) - size
          if (i==2) size.dim1<-1
          else if (i > 2) size.dim1<-dim(p1.freq.all)[1]
          X<-rep(NA,(dif*size.dim1))
          X.mat<-array(X,dim=c(size.dim1,dif))
          if (i==2) {
            p1.freq.all<-c(p1.freq.all,X.mat)
            p2.freq.all<-c(p2.freq.all,X.mat)
          }
          else if (i > 2) {
            p1.freq.all<-cbind(p1.freq.all,X.mat)
            p2.freq.all<-cbind(p2.freq.all,X.mat)
          }	
          p1.freq.all<-rbind(p1.freq.all,p1.freq)
          p2.freq.all<-rbind(p2.freq.all,p2.freq)
        }
        ## if new locus has fewer alleles add NAs to it, if the same number just add it
        else if (length(p1.freq) <= size){
          dif<-size - length(p1.freq)
          p1.freq<-c(p1.freq,rep(NA,dif))
          p2.freq<-c(p2.freq,rep(NA,dif))
          p1.freq.all<-rbind(p1.freq.all,p1.freq)
          p2.freq.all<-rbind(p2.freq.all,p2.freq)
        }
      }
    }
    rownames(p1.freq.all)<-loci.data[,1]
    rownames(p2.freq.all)<-loci.data[,1]
    ## use parental allele frequencies to calculate combos to use
    ## set up data objects
    focal.set<-row.names(parental1)
    combos.touse<-data.frame(
                             obs.delta=numeric(length(focal.set)),
                             max.delta=numeric(length(focal.set)), row.names=focal.set)
    
    maxalleles<-dim(p1.freq.all)[2]
    combos.touse<-data.frame(combos.touse, 
                             combo1=matrix(data=0,nrow=length(focal.set), ncol=maxalleles),
                             combo2=matrix(data=0, nrow=length(focal.set), ncol=maxalleles),
                             spa.c1=numeric(length(focal.set)), 
                             spb.c1=numeric(length(focal.set)),
                             spa.c2=numeric(length(focal.set)), spb.c2=numeric(length(focal.set))
                             )
    ## run test.combinations function to produce combos.touse
    for (i in 1:n.loci){
      combos.touse[i,]<-test.combinations(p1.freq.all[i,],p2.freq.all[i,])
    }
    ## set up array to save allele names; these are needed for indexing in the hybrid index function
    names<-array(dim=c(n.loci,dim(p1.freq.all)[2]))
    rownames(names)<-as.character(loci.data[,1])
    ## this function forces alleles at higher frequency in species two to give a count of 0
    combos.touse<-fixup.combos.touse(combos.touse, loci.data[,1], maxalleles)
    ## calculate count.matrix based on combus.to.use and genotype data
    X<-rep(NA,(n.loci*n.ind))
    count.matrix<-matrix(X,nrow=n.loci,ncol=n.ind)
    for(i in 1:n.loci){
      ## it is important to get the loci names from the parental files for the right loci to be referenced
      a.n <- allele.names(as.genotype(c(parental1[i,],parental2[i,]), allow.partial.missing=TRUE))
      names[i,1:length(a.n)] <- a.n
      ## 16 Oct 09 -- changed the alleles to count to be those from
      ## combo2, which are at high frequency in spa (Parental 1, with
      ## low hybrid index)
      ind.counts<-allele.count(as.genotype(admix.gen[i,],allow.partial.missing=TRUE),
                               names[i,as.numeric(combos.touse[i,(3+maxalleles):(3+2*maxalleles-1)])],
                               na.rm=TRUE)

      ind.counts[is.na(admix.gen[i,])]<-NA
      count.matrix[i,]<-ind.counts      
    }
    rownames(count.matrix)<-loci.data[,1]
    if (is.null(individual.data)==FALSE){
      if (is.na(individual.data[1,2])==FALSE) colnames(count.matrix)<-as.character(individual.data[,2])
    }

    introgress.data<-list(Individual.data=individual.data, Count.matrix=count.matrix,
                          Combos.to.use=combos.touse, Parental1.allele.freq=p1.freq.all,
                          Parental2.allele.freq=p2.freq.all, Alleles=names, Admix.gen=admix.gen)
    return(introgress.data)
  }
}

