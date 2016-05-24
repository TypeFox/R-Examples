#Function: ghap.loadphase
#License: GPLv3 or later
#Modification date: 2 Feb 2016
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Load phased genotypes

ghap.loadphase<-function(
  samples.file,
  markers.file,
  phase.file,
  verbose = TRUE
){
  
  #Check files
  if(file.exists(phase.file) == FALSE){
    stop("Could not find phased genotypes file")
  }
  if(file.exists(markers.file) == FALSE){
    stop("Could not find marker map file")
  }
  if(file.exists(samples.file) == FALSE){
    stop("Could not find sample file")
  }
  
  #Load marker map file
  if(verbose == TRUE){
    cat("\nReading in marker map information... ")
  }
  marker<-read.table(markers.file,header=FALSE,colClasses = c("character","character","numeric","character","character"))
  
  #Check if the map file contains correct dimension
  if(ncol(marker) != 5){
    stop("Marker map contains wrong number of columns (expected 3)")
  }
  
  #Check if the file contains information on a single chromosome
  chr<-unique(marker[,1])
  if(length(chr) != 1){
    stop("Your marker map file contains information on more than one chromosome")
  }
  
  #Check for duplicated marker ids
  if(length(unique(marker[,2])) < nrow(marker)){
    stop("Your marker map file contains duplicated ids!")
  }
  
  #Check if markers are sorted by bp
  if(identical(marker[,3],sort(marker[,3])) == FALSE){
    stop("Markers are not sorted by base pair position")
  }
  
  #Check for duplicated bp
  if(length(unique(marker[,3])) != nrow(marker)){
    warning("Your marker map file contains duplicated ids! Be careful in your analysis!")
  }
  
  #Map passed checks
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nrow(marker), " markers were found for chromosome ",chr,".\n",sep=""))
  }
  
  #Load sample file
  if(verbose == TRUE){
    cat("Reading in sample information... ")
  }
  sample<-read.table(samples.file,header=FALSE,colClasses = "character")
  
  #Check if the sample file contains correct dimension
  if(ncol(sample) != 2){
    stop("Sample file contains wrong number of columns (expected 2)")
  }
  
  #Check for duplicated ids
  if(length(unique(sample[,2])) < nrow(sample)){
    stop("Sample file contains duplicated ids!")
  }
  
  pop <- rep(sample[,1],each=2)
  ids <- rep(sample[,2],each=2)
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nrow(sample), " individuals were found in ", length(unique(pop)), " populations.\n",sep=""))
  }
  
  #Create GHap.phase object
  phase<-NULL
  phase$chr<-chr
  phase$nsamples<-nrow(sample)
  phase$nmarkers<-nrow(marker)
  phase$nsamples.in<-nrow(sample)
  phase$nmarkers.in<-nrow(marker)
  phase$pop<-pop
  phase$id<-ids
  phase$id.in<-rep(TRUE,times=length(phase$id))
  phase$marker<-marker[,2]
  phase$marker.in<-rep(TRUE,times=length(phase$marker))
  phase$bp<-marker[,3]
  phase$A0<-marker[,4]
  phase$A1<-marker[,5]
  if(verbose == TRUE){
    cat("Reading in phased genotypes... (may take a few minutes for large datasets)\n")
  }
  phase$phase<-read.big.matrix(filename = phase.file,sep = " ",header = FALSE,type = "char")
  
  #Check phase file dimensions
  
  if(nrow(phase$phase) != phase$nmarkers & ncol(phase$phase) != 2*phase$nsamples){
    stop("Your phased genotypes file contains wrong dimensions")
  }
  
  #Return ghap object
  class(phase) <- "GHap.phase"
  if(verbose == TRUE){
    cat("Your GHap.phase object was successfully loaded without apparent errors.\n\n")
  }
  return(phase)
  
}