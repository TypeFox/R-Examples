#Function: ghap.subsethaplo
#License: GPLv3 or later
#Modification date: 2 Feb 2016
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Subset a GHap.haplo object

ghap.subsethaplo<-function(
  haplo,
  ids,
  alleles,
  verbose = TRUE
){
  
  #Check if haplo is a GHap.haplo object
  if(class(haplo) != "GHap.haplo"){
    stop("Argument haplo must be a GHap.haplo object.")
  }
  
  #Log message
  if(verbose == TRUE){
    cat("\nSubsetting ", length(ids), " individuals and ", length(which(alleles)), " haplotype alleles... ",sep="")
  }
  
  #Subsetting
  haplo$allele.in<-alleles
  haplo$id.in<-haplo$id %in% ids
  haplo$nsamples.in<-length(which(haplo$id.in))
  haplo$nalleles.in<-length(which(haplo$allele.in))
  
  #Log message
  if(verbose == TRUE){
    cat("Done.\n")
    cat("Final data contains ", haplo$nsamples.in, " individuals and ", haplo$nalleles.in, " haplotype alleles.\n\n",sep="")
  }
  
  #Output object
  return(haplo)
  
}