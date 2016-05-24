#Function: ghap.outhaplo
#License: GPLv3 or later
#Modification date: 2 Feb 2016
#Written by: Marco Milanesi & Yuri Tani Utsunomiya
#Contact: marco.milanesi.mm@gmail.com, ytutsunomiya@gmail.com
#Description: Output a GHap.haplo object as a file

ghap.outhaplo<-function(haplo,outfile,only.active.markers=TRUE,only.active.samples=TRUE,verbose=TRUE){
  
  #Check if inactive markers and samples should be reactived
  if(only.active.markers == FALSE){
    haplo$allele.in <- rep(TRUE,times=haplo$nalleles)
  }
  if(only.active.samples == FALSE){
    haplo$id.in <- rep(TRUE,times=haplo$nsamples)
  }
  
  #Insert suffix to outfile
  hapgenotypes <- paste(outfile,"hapgenotypes",sep=".")
  hapalleles <- paste(outfile,"hapalleles",sep=".")
  hapsamples <- paste(outfile,"hapsamples",sep=".")
  
  #Check if output will overwrite existing files before opening connection
  if(file.exists(hapsamples) == TRUE){
    stop(paste("File", hapsamples, "already exists!"))
  }
  if(file.exists(hapgenotypes) == TRUE){
    stop(paste("File", hapgenotypes, "already exists!"))
  }
  if(file.exists(hapalleles) == TRUE){
    stop(paste("File", hapalleles, "already exists!"))
  } 
  
  
  #Prepare hapsample file
  if(verbose == TRUE){
    cat("\nPreparing ", hapsamples, "... ", sep="")
  }
  id <- haplo$id[haplo$id.in]
  pop <- haplo$pop[haplo$id.in]
  
  #Output hapsamples file
  write.table(x = cbind(pop,id),file = hapsamples,quote = FALSE,row.names = FALSE,col.names=FALSE)
  if(verbose == TRUE){
    cat("Done.\n")
  }
  
  
  #Prepare hapalleles file
  if(verbose == TRUE){
    cat("Preparing ", hapalleles, "... ", sep="")
  }
  map <- NULL
  map$BLOCK <- haplo$block[haplo$allele.in]
  map$CHR <- haplo$chr[haplo$allele.in]
  map$BP1 <- haplo$bp1[haplo$allele.in]
  map$BP2 <- haplo$bp2[haplo$allele.in]
  map$ALLELE <- haplo$allele[haplo$allele.in]
  map <- as.data.frame(map)
  
  #Output hapalleles file
  write.table(x = map, file = hapalleles, quote = FALSE,sep=" ", row.names = FALSE, col.names=FALSE)
  if(verbose == TRUE){
    cat("Done.\n")
  }
  
  #Prepare hapgenotypes file
  if(verbose == TRUE){
    cat("Preparing ", hapgenotypes, "... ", sep="")
  }
  write.big.matrix(x = as.big.matrix(haplo$genotypes[haplo$allele.in,haplo$id.in]), filename = hapgenotypes, sep = " ", row.names = FALSE, col.names = FALSE)
  if(verbose == TRUE){
    cat("Done.\n\n", sep="")   
  }
  
}

