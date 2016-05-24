

`genomesim` <-
function(pedigree, founders=NULL, positions=NULL, initHe=NULL, 
            mutationType=NULL, mutationRate=NULL, phenotyped=NULL, 
            founderHaplotypes=NULL, genotyped=NULL, returnG='n',initFreqs=NULL){


original.order<-order(positions)
positions<-positions[original.order]
initHe<-initHe[original.order]
mutationType<-mutationType[original.order]
mutationRate<-mutationRate[original.order]
founderHaplotypes<-founderHaplotypes[,original.order]


if(all.equal(original.order,1:length(original.order))==FALSE&is.null(initFreqs)==FALSE) {
	stop("ordering of initFreqs will fail")
}


"extractA"<-function(G,marker.type="MSW"){
## stollen from Jarrod to avoid
## dependency on MasterBayes

   if(is.genotype(G[[1]])==FALSE & is.genotypeD(G[[1]])==FALSE){
     if("id"%in%colnames(G)){
       G<-G[,-which(colnames(G)=="id")]
     }
     if("categories"%in%colnames(G)){
       G<-G[,-which(colnames(G)=="categories")]
     }
     G<-genotype.list(G,marker.type=marker.type)
   }

  A<-lapply(G, function(x){summary(x)$allele.freq[,"Proportion"][1:nallele(x)]})
  if(marker.type=="AFLP" | is.genotypeD(G[[1]])){
  A<-lapply(A, function(x){
       if(any(x==1)){
         warning("some loci are monomorphic: setting allele frequencies to 0.9999 and 0.0001")
         x[which(x==1)]<-0.9999
         x[which(x==0)]<-0.0001
       }else{
         x<-x
       }
       x
     })
  }
  A
}


  for(x in 1:3) pedigree[,x]<-as.character(pedigree[,x])

  if(any(match(pedigree[,2], pedigree[,1])>match(pedigree[,1], pedigree[,1]), na.rm=T))
                 stop("Dams appearing after their offspring: use fixPedigree")
  if(any(match(pedigree[,3], pedigree[,1])>match(pedigree[,1], pedigree[,1]), na.rm=T))
                 stop("Sires appearing after their offspring: use fixPedigree")

  if(is.null(phenotyped)==FALSE&&is.numeric(phenotyped)&&length(phenotyped)>length(pedigree[,1]))
                 stop("Length of phenotype availability indicator vector incompatible with length of pedigree")

  ## assumes that 'true pedigree has been supplied' and that conditions
  ## of fixPedigreep()

  ## haplotypes can be supplied for founders, in which case they will be assigned in the order provided
  ## i.e., at least twice as many haplotypes as founders are required.
  founderHap<-1

  ## mutationType can be 'cIAM','dIAM','Micro','Dom'
  mutation<-function(currentAllele,mutationType,mu){
    mutAllele<-currentAllele
    if(rbinom(1,1,mu)==1){
      if(mutationType=='cIAM'|mutationType=='dIAM') mutAllele<-rnorm(1,0,1)
      if(mutationType=='Micro'){
        if(rbinom(1,1,0.5)==1) {mutAllele<-currentAllele-1} else {mutAllele<-currentAllele+1}
      }
      if(mutationType=='dom') mutAllele<-abs(currentAllele-1)
    }
    mutAllele
  }

  ## while the continuous IAM model will work well for mutations, it will cause
  ## A = N.founders if applied directly.  So I'll make a matrix of loci times
  ## times 1000, fill and fill it with founder alleles are appropriate frequencies
  if(is.null(initFreqs)&is.null(founderHaplotypes)){
    founderIAMalleles<-matrix(NA,length(positions),1000)
    for(x in 1:length(positions)){
      A<-as.integer(1/(1-initHe[x])+0.9999)
      A.index<-as.integer(runif(1000,0,1/(1-initHe[x])))+1
      A.continuous.values<-rnorm(A,0,1)
      founderIAMalleles[x,]<-A.continuous.values[A.index]
    }
  }

  ## mutationType can be 'cIAM','dIAM','Micro','Dom'
  rndAllele<-function(loc){
    allele<-NULL
    if(is.null(initFreqs)) {
      if(mutationType[loc]=='cIAM'|mutationType[loc]=='dIAM') allele<-founderIAMalleles[loc,as.integer(runif(1,0,1000))+1]
      if(mutationType[loc]=='Micro'){
        effectiveAlleles<-1/(1-initHe[loc])
        allele<-as.integer(runif(1,0,effectiveAlleles))+200
      }
      if(mutationType[loc]=='Dom'){
        effectiveAlleles<-1/(1-initHe[loc])
        allele<-as.integer(runif(1,0,effectiveAlleles))
      }
    }else{
      i<-runif(1,0,1)
      a<-0
      while(i>0) { a<-a+1; i<-i-initFreqs[[loc]][a]; }
      allele<-as.numeric(names(initFreqs[[loc]][a]))
      if(mutationType[loc]=='cIAM') allele<-as.numeric(names(initFreqs[[loc]][a]))
    }
    allele
  }

  ## to get the right bit go to (individual-1)*2+chromosome
  genomes<-matrix(NA,length(pedigree[,1])*2,length(positions))
  
  Kosambi.c<-function(m){
  	f<-function(m,c){m-0.25*log((1+2*c)/(1-2*c))}
  	c<-m
  	if(m>0) c<- uniroot(f,c(0,0.5),m=m)$root
  	c
  }

  segregate<-function(parent){
  	# map distances between adjascent loci
  	map.intervals<-positions[2:length(positions)]-
  	                  positions[1:(length(positions)-1)]
  	# Haldane-based expected crossover freq 
#  	obs.cross.frac<-((1-exp(-2*map.intervals/100))/2)
    # Kosambi's mapping function
    obs.cross.frac<-apply(X=as.matrix(map.intervals)/100,MARGIN=1,FUN=Kosambi.c)

  	# observed recombinations
  	recomb<-rbinom(length(obs.cross.frac),1,obs.cross.frac)

    p<-genomes[((parent-1)*2+1):((parent-1)*2+2),]
    
    offspringHaplotype<-array(dim=length(positions))
    chromosome<-sample(1:2,1)
    offspringHaplotype[1]<-p[chromosome,1]
    for(x in 1:length(map.intervals)){
    	if(recomb[x]==1&chromosome==1){ chromosome<-2
    		}else{
    	if(recomb[x]==1&chromosome==2) chromosome<-1
    	}
    	offspringHaplotype[x+1]<-p[chromosome,x+1]
    }
  	
    offspringHaplotype
  }

  cat(paste("Processing pedigree...",'\n')); flush.console();
  cat(paste("0%                     50%                     100%",'\n')); flush.console();
  cat(paste("|                       |                       |",'\n')); flush.console();
  mark.at<-(length(pedigree[,1])-1)/50
  count<-0

  ## now to go through all individuals
  for(x in 1:length(pedigree[,1])){
    count<-count+1
    if(count>mark.at){ count<-0; cat("-"); flush.console(); }

    ## assign HW-assumed genotypes if founder
#    if(founders[x]==1){
#      genomes[(x-1)*2+1,]<-mapply(rndAllele,loc=1:length(positions))
#      genomes[(x-1)*2+2,]<-mapply(rndAllele,loc=1:length(positions))

#    } else {
    ## ... otherwise base genotypes on parents

      if(is.na(pedigree$dam[x])) {
        if(is.null(founderHaplotypes)) { genomes[(x-1)*2+1,]<-mapply(rndAllele,loc=1:length(positions))
        } else {
          genomes[(x-1)*2+1,]<-founderHaplotypes[founderHap,]
          founderHap<-founderHap+1
        }
      } else {
          mum<-which(as.character(pedigree$id)==as.character(pedigree$dam[x]),arr.ind=TRUE)
      ## segregation
        genomes[(x-1)*2+1,]<-segregate(mum)
        ## mutation
        genomes[(x-1)*2+1,]<-mapply(mutation,currentAllele=genomes[(x-1)*2+1,],mutationType=mutationType,mu=mutationRate)
      }
      if(is.na(pedigree$sire[x])) {
        if(is.null(founderHaplotypes)) { genomes[(x-1)*2+2,]<-mapply(rndAllele,loc=1:length(positions))
        } else {
          genomes[(x-1)*2+2,]<-founderHaplotypes[founderHap,]
          founderHap<-founderHap+1
        }
      } else {
        dad<-which(as.character(pedigree$id)==as.character(pedigree$sire[x]),arr.ind=TRUE)
        ## segregation
        genomes[(x-1)*2+2,]<-segregate(dad)
        ## mutation
        genomes[(x-1)*2+2,]<-mapply(mutation,currentAllele=genomes[(x-1)*2+2,],mutationType=mutationType,mu=mutationRate)
      }

  }
  cat(paste('\n',"...done.",'\n')); flush.console();

  cat(paste('\n',"Calculating phenotypes...")); flush.console();
  allelicEffect<-1

  genes<-genomes[subset(1:length(genomes[,1]),(1:length(genomes[,1]))%%2==0),subset(1:length(positions),mutationType=='cIAM')]+
              genomes[subset(1:length(genomes[,1]),(1:length(genomes[,1]))%%2==1),subset(1:length(positions),mutationType=='cIAM')]

  if(is.null(dim(genes)[2])==FALSE&&dim(genes)[2]>1) {
     phenotypes<-as.data.frame(cbind(pedigree$id,rowSums(genes)*(1*allelicEffect)))
  } else {
     phenotypes<-as.data.frame(cbind(pedigree$id,genes*(1*allelicEffect)))
  }
  if(is.character(phenotyped)) phenotypes<-phenotypes[which(phenotypes[,1]%in%phenotyped),]
  if(is.numeric(phenotyped)) phenotypes<-subset(phenotypes,phenotyped>0)
  cat(paste("done.",'\n')); flush.console();

  markerData<-NULL

  if(length(subset(mutationType,mutationType!='cIAM'))>0){
    cat(paste('\n',"Tabulating marker genotypes...")); flush.console();
    for(x in 1:length(positions)) {if (mutationType[x]=='dIAM') genomes[,x]<-as.numeric(as.factor(genomes[,x]))}
    markers<-genomes[,subset(1:length(positions),mutationType!='cIAM')]
    doubleids<-as.vector(as.matrix(as.data.frame(rbind(pedigree$id,pedigree$id))))
    markers<-as.data.frame(cbind(doubleids,markers))
    markers<-markers[which(markers[,1]%in%genotyped),]
    markerData<-markers
    idsformarkers<-markers[,1]
    markers<-as.matrix(markers[,2:length(markers[1,])])
    cat(paste("done.",'\n')); flush.console();

    mt<-subset(mutationType, mutationType!='cIAM')
    for(x in 1:length(mt)) {
      if(length(table(markers[,x]))>2) { mt[x]<-'M' } else { mt[x]<-'S' }
    }
  }


#  if(FALSE){
#    write(length(markers[,1])/2, "~/Desktop/textGenomeSim.out", append=TRUE)
#    write(length(markers[1,]), "~/Desktop/textGenomeSim.out", append=TRUE)
#    markerPositions<-subset(positions, mutationType!='cIAM')
#    p<-array(dim=length(markerPositions)+1); p[2:length(p)]<-markerPositions; p[1]<-"P";
#    write(p, "~/Desktop/textGenomeSim.out", ncolumns = length(p), append=TRUE)
#    write(mt, "~/Desktop/textGenomeSim.out", ncolumns = length(mt), append=TRUE, sep="")
#    for(x in 1:(length(markers[,1])/2)){
#      write(idsformarkers[(x-1)*2+1], "~/Desktop/textGenomeSim.out", append=TRUE)
#      write(markers[(x-1)*2+1,], "~/Desktop/textGenomeSim.out", ncolumns = length(markers[,1]), append=TRUE)
#      write(markers[(x-1)*2+2,], "~/Desktop/textGenomeSim.out", ncolumns = length(markers[,1]), append=TRUE)
#    }
#  }

## calculate pairwise similarities a la Hayes and Goddard
#  S<-matrix(NA,(length(genes[,1])/2),(length(genes[,1])/2))

#  for(x in 1:(length(S[1,])-1)){
#cat(paste(x,'\n')); flush.console();
#    for(y in x:length(S[,1])){
#      sim=0;
#      for(g in 1:length(genes[1,])){
#        obs<-length(table(c(genes[((x-1)*2+1):((x-1)*2+2),g],genes[((y-1)*2+1):((y-1)*2+2),g])))
#        if(obs<4) sim<-sim+1/(2^(obs-1))
#      }
#      S[x,y]=sim/length(genes[1,])
#    }
#  }
#  S<-(S-min(S))/1-min(S)
#  output<-list(Phenotypes=phenotypes,MarkerData=markerData,MarkerRelatedness=S)

  output<-list(Phenotypes=phenotypes,MarkerData=markerData)
  if(returnG=='y') output<-list(Phenotypes=phenotypes,MarkerData=markerData,genomes=as.data.frame(cbind(doubleids,genomes)))

  output
}











