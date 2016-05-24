`microsim` <-
function(pedigree,genFreqs=NULL,genotypesSample=NULL,knownGenotypes=NULL,records=NULL,eRate1=0,eRate2=0,eRate3=0){

  if(is.null(genotypesSample)==FALSE&&(names(genotypesSample)[1]=="id"|names(genotypesSample)[1]=="ID")) genotypesSample<-genotypesSample[,-1]
  if(is.null(genotypesSample)==FALSE&is.null(genFreqs)==TRUE) genFreqs<-extractA(genotypesSample)


  loci<-length(names(genFreqs))

if(is.null(records)==FALSE&&is.data.frame(records)==FALSE&&is.matrix(records)==FALSE&&length(records)!=loci&&length(records)!=1)  
           stop("Dimension of records incompatible with number of loci")
if(is.null(records)==FALSE&&(is.data.frame(records)|is.matrix(records))&&dim(records)[1]!=length(pedigree[,1]))  
           stop("Dimension of records incompatible with the size of the pedigree")
if(is.null(records)==FALSE&&(is.data.frame(records)|is.matrix(records))&&dim(records)[2]!=loci)  
           stop("Dimension of records incompatible with number of loci")

  if(length(eRate1)!=1&eRate1!=loci) stop("Length of eRate1 incompatible with number of loci")
  if(length(eRate2)!=1&eRate2!=loci) stop("Length of eRate2 incompatible with number of loci")
  if(length(eRate3)!=1&eRate3!=loci) stop("Length of eRate3 incompatible with number of loci")

  microGen<-as.data.frame(matrix(NA,length(pedigree[,1]),loci*2+1))
  microGen[,1]<-as.character(pedigree[,1])
  row.names(microGen)<-as.character(pedigree[,1])
  microGen$genotyped<-0

  randAllele<-function(loc,freqs){
    i<-runif(1,0,1)
    a<-0
    while(i>0) { a<-a+1; i<-i-freqs[[loc]][a]; }
    names(freqs[[loc]][a])
  }


  assignGenotype<-function(id,mum,dad,f){
    maternalGenotype<-NULL
    paternalGenotype<-NULL

    if(is.na(mum)==TRUE) {
      maternalGenotype<-array(dim=loci*2)
      for(l in 1:loci) { maternalGenotype[(l-1)*2+1]<-randAllele(l,f); maternalGenotype[(l-1)*2+2]<-randAllele(l,f); }
    } else {
      if(microGen[mum,dim(microGen)[2]]==1)
               maternalGenotype<-subset(microGen,microGen[,1]==mum)[1,2:(loci*2+1)]
    }

    if(is.na(dad)==TRUE) {
      paternalGenotype<-array(dim=loci*2)
      for(l in 1:loci) { paternalGenotype[(l-1)*2+1]<-randAllele(l,f); paternalGenotype[(l-1)*2+2]<-randAllele(l,f); }
    } else {
      if(microGen[dad,dim(microGen)[2]]==1)
               paternalGenotype<-subset(microGen,microGen[,1]==dad)[1,2:(loci*2+1)]
    }

    if(is.null(maternalGenotype)==FALSE&is.null(paternalGenotype)==FALSE) {
      microGen[id,dim(microGen)[2]]<<-1
      for(l in 1:loci){
        maternalAllele<-rbinom(1,1,0.5)
        paternalAllele<-rbinom(1,1,0.5)
        microGen[id,(l-1)*2+2]<<-maternalGenotype[(l-1)*2+1+maternalAllele]
        microGen[id,(l-1)*2+3]<<-paternalGenotype[(l-1)*2+1+paternalAllele]
      }
    }
  }

  while(prod(microGen$genotyped)==0){
    for(x in 1:dim(microGen)[1]) 
         assignGenotype(as.character(pedigree$id[x]),as.character(pedigree$dam[x]),as.character(pedigree$sire[x]),genFreqs)
  }

  trueMicroGenotypes<-microGen[-c(1,dim(microGen)[2])]
  for(l in 1:loci) names(trueMicroGenotypes)[(l-1)*2+1]<-names(genFreqs[l])
  for(l in 1:loci) names(trueMicroGenotypes)[(l-1)*2+2]<-paste(names(genFreqs[l]),"_b",sep="")


  observedMicroGenotypes<-trueMicroGenotypes

  # make some genotypic errors

  # genotype substitutions
  if(sum(eRate1)>0) {
    if(length(eRate1)==1){
      for(i in 1:dim(observedMicroGenotypes)[1]){
        for(l in 1:loci) {
          if(runif(1,0,1)<eRate1) {
            observedMicroGenotypes[i,(l-1)*2+1]<-randAllele(l,genFreqs)
            observedMicroGenotypes[i,(l-1)*2+2]<-randAllele(l,genFreqs)
          }
        }
      }
    }
    if(length(eRate1)==loci){
      for(i in 1:dim(observedMicroGenotypes)[1]){
        for(l in 1:loci) {
          if(runif(1,0,1)<eRate1[l]) {
            observedMicroGenotypes[i,(l-1)*2+1]<-randAllele(l,genFreqs)
            observedMicroGenotypes[i,(l-1)*2+2]<-randAllele(l,genFreqs)
          }
        }
      }
    }
  }
  # allele substitutions
  if(sum(eRate2)>0) {
    if(length(eRate2)==1){
      for(i in 1:dim(observedMicroGenotypes)[1]){
        for(l in 1:loci) {
          if(runif(1,0,1)<eRate2) observedMicroGenotypes[i,(l-1)*2+1]<-randAllele(l,genFreqs)
          if(runif(1,0,1)<eRate2) observedMicroGenotypes[i,(l-1)*2+2]<-randAllele(l,genFreqs)
        }
      }
    }
    if(length(eRate2)==loci){
      for(i in 1:dim(observedMicroGenotypes)[1]){
        for(l in 1:loci) {
          if(runif(1,0,1)<eRate2[l]) observedMicroGenotypes[i,(l-1)*2+1]<-randAllele(l,genFreqs)
          if(runif(1,0,1)<eRate2[l]) observedMicroGenotypes[i,(l-1)*2+2]<-randAllele(l,genFreqs)
        }
      }
    }
  }
  # large allele dropouts
  if(sum(eRate3)>0) {
    if(length(eRate3)==1){
      for(i in 1:dim(observedMicroGenotypes)[1]){
        for(l in 1:loci) {
          if(runif(1,0,1)<eRate3) {
            smallerAllele<-min(c(trueMicroGenotypes[i,(l-1)*2+1],trueMicroGenotypes[i,(l-1)*2+2]))
            observedMicroGenotypes[i,(l-1)*2+1]<-smallerAllele
            observedMicroGenotypes[i,(l-1)*2+2]<-smallerAllele
          }
        }
      }
    }
    if(length(eRate3)==loci){
      for(i in 1:dim(observedMicroGenotypes)[1]){
        for(l in 1:loci) {
          if(runif(1,0,1)<eRate3[l]) {
            smallerAllele<-min(c(trueMicroGenotypes[i,(l-1)*2+1],trueMicroGenotypes[i,(l-1)*2+2]))
            observedMicroGenotypes[i,(l-1)*2+1]<-smallerAllele
            observedMicroGenotypes[i,(l-1)*2+2]<-smallerAllele
          }
        }
      }
    }
  }

  # cut out some missing records
  genAvailIndex<-NULL
  if(is.null(records)==FALSE) {
    recProbs<-matrix(0,dim(trueMicroGenotypes)[1],dim(trueMicroGenotypes)[2])
    if(dim(as.data.frame(records))[1]==1&dim(as.data.frame(records))[2]==1) 
                          recProbs<-matrix(records,dim(trueMicroGenotypes)[1],dim(trueMicroGenotypes)[2])
    if(length(records)==loci) for(l in 1:loci) recProbs[,(l-1)*2+1]<-records[l]
    if(dim(as.data.frame(records))[1]==length(pedigree[,1])&dim(as.data.frame(records))[2]==loci) 
                         for(l in 1:loci) recProbs[,(l-1)*2+1]<-records[,l]
    rbinfunc<-function(x) (rbinom(1,1,x))
    genAvailIndex<-apply(recProbs,c(1,2),rbinfunc)
    for(l in 1:loci) genAvailIndex[,(l-1)*2+2]<-genAvailIndex[,(l-1)*2+1]
  }

  if(is.null(genAvailIndex)==FALSE) {
    missingtozeros<-matrix(as.numeric(as.matrix(observedMicroGenotypes)),dim(observedMicroGenotypes)[1],dim(observedMicroGenotypes)[2],byrow=FALSE)*genAvailIndex
    zeromissingfunc<-function(x) {res<-x; if(x==0) res<-NA; res;}
    withMissing<-as.data.frame(apply(missingtozeros,c(1,2),zeromissingfunc))
    names(withMissing)<-names(observedMicroGenotypes)
    row.names(withMissing)<-row.names(observedMicroGenotypes)
    observedMicroGenotypes<-withMissing
  }

  list(trueGenotypes=trueMicroGenotypes, observedGenotypes=observedMicroGenotypes)
}

