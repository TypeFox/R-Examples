#Function: ghap.loadhaplo
#License: GPLv3 or later
#Modification date: 2 Feb 2016
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Load haplotype genotypes

ghap.loadhaplo<-function(
  hapsamples.file,
  hapalleles.file,
  hapgenotypes.file,
  verbose = TRUE
){
  
  #Check files
  if(file.exists(hapgenotypes.file) == FALSE){
    stop("Could not find haplotype genotypes file")
  }
  if(file.exists(hapalleles.file) == FALSE){
    stop("Could not find haplotype alleles file")
  }
  if(file.exists(hapsamples.file) == FALSE){
    stop("Could not find sample file")
  }
  
  #Load haplotype alleles file
  if(verbose == TRUE){
    cat("\nReading in haplotype allele information... ")
  }
  hapalleles<-read.table(hapalleles.file,header=FALSE,colClasses = c("character","character","numeric","numeric","character"))
  
  #Check if the haplotype alleles file contains correct dimension
  if(ncol(hapalleles) != 5){
    stop("Haplotype alleles file contains wrong number of columns (expected 5)")
  }else{
    if(verbose == TRUE){
      cat("Done.\n")
      cat(paste("A total of ", nrow(hapalleles), " haplotype alleles were found.\n",sep=""))
    }
  }
  
  #Load sample file
  if(verbose == TRUE){
    cat("Reading in sample information... ")
  }
  sample<-read.table(hapsamples.file,header=FALSE,colClasses = "character")
  
  #Check if the sample file contains correct dimension
  if(ncol(sample) != 2){
    stop("Sample file contains wrong number of columns (expected 2)")
  }
  
  #Check for duplicated ids
  if(length(unique(sample[,2])) < nrow(sample)){
    stop("Sample file contains duplicated ids!")
  }
  
  pop<-sample[,1]
  ids<-sample[,2]
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nrow(sample), " individuals were found in ", length(unique(pop)), " populations.\n",sep=""))
  }
  
  #Create GHap.hap object
  hap<-NULL
  hap$nsamples<-nrow(sample)
  hap$nalleles<-nrow(hapalleles)
  hap$nsamples.in<-nrow(sample)
  hap$nalleles.in<-nrow(hapalleles)
  hap$pop<-pop
  hap$id<-ids
  hap$id.in<-rep(TRUE,times=length(hap$id))
  hap$chr<-hapalleles[,2]
  hap$block<-hapalleles[,1]
  hap$bp1<-hapalleles[,3]
  hap$bp2<-hapalleles[,4]
  hap$allele<-hapalleles[,5]
  hap$allele.in<-rep(TRUE,times=length(hap$allele))
  if(verbose == TRUE){
    cat("Reading in haplotype genotypes... (may take a few minutes for large datasets)\n")
  }
  hap$genotypes<-read.big.matrix(filename = hapgenotypes.file,sep = " ",header = FALSE,type = "char")
  
  #Check haplotype genotypes file dimensions
  if(nrow(hap$genotypes) != hap$nalleles & ncol(hap$genotypes) != 2*hap$nsamples){
    stop("Your haplotype genotypes file contains wrong dimensions")
  }
  
  #Return GHap.haplo object
  class(hap) <- "GHap.haplo"
  if(verbose == TRUE){
    cat("Your GHap.haplo object was successfully loaded without apparent errors.\n\n")
  }
  return(hap)
  
}
