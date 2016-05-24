#Function: ghap.mergephase
#License: GPLv3 or later
#Modification date: 13 Apr 2016
#Written by: Marco Milanesi & Yuri Tani Utsunomiya
#Contact: marco.milanesi.mm@gmail.com, ytutsunomiya@gmail.com
#Description: Merge GHap.phase objects 

ghap.mergephase<-function(
  phase.1,  phase.2,  only.active.markers=TRUE, only.active.samples=TRUE, verbose=TRUE
){
  
  phase1 <- phase.1
  phase2 <- phase.2
  rm(phase.1,phase.2)
  
  #Check if phase1 is a GHap.phase object
  if(class(phase1) != "GHap.phase"){
    stop("Argument phase1 must be a GHap.phase object.")
  }
  #Check if phase2 is a GHap.phase object
  if(class(phase2) != "GHap.phase"){
    stop("Argument phase2 must be a GHap.phase object.")
  }
  
  #Check if inactive markers and samples should be reactived
  if(only.active.markers == FALSE){
    phase1$marker.in <- rep(TRUE,times=phase1$nmarkers)
    phase2$marker.in <- rep(TRUE,times=phase2$nmarkers)
    phase1$nmarkers.in <- phase1$nmarkers
    phase2$nmarkers.in <- phase2$nmarkers
  }
  if(only.active.samples == FALSE){
    phase1$id.in <- rep(TRUE,times=2*phase1$nsamples)
    phase2$id.in <- rep(TRUE,times=2*phase2$nsamples)
    phase1$nsamples.in <- phase1$nsamples
    phase2$nsamples.in <- phase2$nsamples
  }
  
  #Check if markers are the same and have the same order in the two objects
  if(identical(phase1$marker[phase1$marker.in],phase2$marker[phase2$marker.in]) == FALSE | (phase1$chr==phase2$chr) == FALSE){
    stop("Markers in the two GHap.phase objects are not the same or are not in the same order.")
  }
  
  #Check if there is some individual duplicate
  if(length(which(phase1$id[phase1$id.in] %in% phase2$id[phase2$id.in])) != 0){
    stop("Duplicated IDs are not allowed.")
  }
  
  
  #Create the merged GHap.phase object
  if(verbose == TRUE){
    cat("\nCreating the new GHap.phase object... ")
  }
  phase<-NULL
  phase$chr<-phase1$chr
  phase$nsamples<-(phase1$nsamples.in+phase2$nsamples.in)
  phase$nmarkers<-phase1$nmarkers.in
  phase$nsamples.in<-(phase1$nsamples.in+phase2$nsamples.in)
  phase$nmarkers.in<-phase1$nmarkers.in
  phase$pop<-c(phase1$pop[phase1$id.in],phase2$po[phase2$id.in])
  phase$id<-c(phase1$id[phase1$id.in],phase2$id[phase2$id.in])
  phase$id.in<-rep(TRUE,(2*phase$nsamples.in))
  phase$marker<-phase1$marker[phase1$marker.in]
  phase$marker.in<-rep(TRUE,phase$nmarkers.in)
  phase$bp<-phase1$bp[phase1$marker.in]
  phase$A0<-phase1$A0[phase1$marker.in]
  phase$A1<-phase1$A1[phase1$marker.in]
  tmp1 <- as.matrix(phase1$phase[phase1$marker.in,phase1$id.in]) # export to a matrix first phase data
  tmp2 <- as.matrix(phase2$phase[phase2$marker.in,phase2$id.in]) # export to a matrix second phase data
  options(bigmemory.typecast.warning=FALSE)
  phase$phase <- as.big.matrix(cbind(tmp1,tmp2), type="char")
  rm(tmp1,tmp2)
  
  #Return ghap object
  class(phase) <- "GHap.phase"
  if(verbose == TRUE){
    cat("Done.\n")
    cat("Your GHap.phase object was successfully merged without apparent errors.\n\n")
  }
  return(phase)
}
