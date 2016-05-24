#Function: ghap.subsetphase
#License: GPLv3 or later
#Modification date: 2 Feb 2016
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Subset a GHap.phase object

ghap.subsetphase<-function(
  phase,
  ids,
  markers,
  verbose = TRUE
){
  
  #Check if phase is a GHap.phase object
  if(class(phase) != "GHap.phase"){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  #Log message
  if(verbose == TRUE){
    cat("\nSubsetting ", length(ids), " individuals and ", length(markers), " markers... ",sep="")
  }
  
  #Subsetting
  phase$marker.in<-phase$marker %in% markers
  phase$id.in<-phase$id %in% ids
  phase$nsamples.in<-length(which(phase$id.in))/2
  phase$nmarkers.in<-length(which(phase$marker.in))
  
  #Log message
  if(verbose == TRUE){
    cat("Done.\n")
    cat("Final data contains ", phase$nsamples.in, " individuals and ", phase$nmarkers.in, " markers.\n ",sep="")
  }
  
  #Output object
  return(phase)
  
}