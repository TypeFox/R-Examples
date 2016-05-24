test.genotypes <-
function(admix.gen=NULL,loci.data=NULL,parental1=NULL,
                         parental2=NULL){
  ##the genetics library is needed for this function
  ##library(genetics)

  ## set up files test for admixed alleles not sampled from
  ## the parental species
  dropped.data<-array(rep(0,dim(loci.data)[1]),dim=c(dim(loci.data)[1],1))
  rownames(dropped.data)<-loci.data[,1]
  n.loci<-dim(admix.gen)[1]
  n.ind<-dim(admix.gen)[2]
  for (i in 1:n.loci){
    parental.alleles<-allele.names(as.genotype(c(parental1[i,],parental2[i,]),
                                               allow.partial.missing=TRUE))
    for (j in 1:n.ind){
      if (is.na(admix.gen[i,j])==FALSE){
        present<-1
        admix.alleles<-allele.names(as.genotype(admix.gen[i,j],allow.partial.missing=TRUE))
        for (k in 1:length(admix.alleles)){
          if (sum(admix.alleles[k]==parental.alleles)==0) present<-0
        }
        if (present==0) {
          admix.gen[i,j]<-NA
          dropped.data[i,1]<-dropped.data[i,1]+1
        }
      }
    } 
  }
  ## only print dropped.data if some data were dropped
  if (sum(dropped.data[,1]) > 0){
    cat("Genotype data for the following number of individuals were dropped", 
         "because admixed individuals had alleles not encountered in either",
         "parental population", fill=TRUE)
    print (dropped.data)
  }  
  return (admix.gen)
}

