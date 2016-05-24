#Function: ghap.out
#License: GPLv3 or later
#Modification date: 2 Feb 2016
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Output a GHap.phase object as a file

ghap.outphase<-function(
  phase,
  outfile,
  only.active.markers = TRUE,
  only.active.samples = TRUE,
  verbose = TRUE
){
  
  #Check if inactive markers and samples should be reactived
  if(only.active.markers == FALSE){
    phase$marker.in <- rep(TRUE,times=phase$nmarkers)
  }
  if(only.active.samples == FALSE){
    phase$id.in <- rep(TRUE,times=2*phase$nsamples)
  }
  
  #Check if output will ovewrite existing files
  samplefile <- paste(outfile,"samples",sep=".")
  mapfile <- paste(outfile,"markers",sep=".")
  phasefile <- paste(outfile,"phase",sep=".")
  
  if(file.exists(samplefile) == TRUE){
    stop("The samples file already exists!")
  }
  if(file.exists(mapfile) == TRUE){
    stop("The markers file already exists!")
  }
  if(file.exists(phasefile) == TRUE){
    stop("The phased genotypes file already exists!")
  }
  
  #Prepare .markers file
  if(verbose == TRUE){
    cat("\nPreparing ", outfile, ".markers... ", sep="")
  }
  map <- NULL
  map$CHR <- phase$chr
  map$SNP <- phase$marker[phase$marker.in]
  map$BP <- phase$bp[phase$marker.in]
  map$A0 <- phase$A0[phase$marker.in]
  map$A1 <- phase$A1[phase$marker.in]
  map <- as.data.frame(map)
  write.table(x = map, file = mapfile,quote = FALSE,sep=" ", row.names = FALSE, col.names=FALSE)
  if(verbose == TRUE){
    cat("Done.\n")
  }
  
  #Prepare .samples file
  if(verbose == TRUE){
    cat("Preparing ", outfile, ".samples... ", sep="")
  }
  samp <- NULL
  ids <- phase$id[phase$id.in]
  ids <- ids[1:length(ids) %% 2 == 0]
  pops <- phase$pop[phase$id.in]
  pops <- pops[1:length(pops) %% 2 == 0]
  samp$POP <- pops
  samp$ID <- ids
  samp <- as.data.frame(samp)
  write.table(x = samp, file = samplefile,quote = FALSE,sep=" ", row.names = FALSE, col.names=FALSE)
  if(verbose == TRUE){
    cat("Done.\n")
  }
  
  #Prepare .phase file
  if(verbose == TRUE){
    cat("Preparing ", outfile, ".phase... ", sep="")   
  }
  write.big.matrix(x = as.big.matrix(phase$phase[phase$marker.in,phase$id.in]), filename = phasefile, sep = " ", row.names = FALSE, col.names = FALSE)
  if(verbose == TRUE){
    cat("Done.\n\n", sep="")   
  }
}

