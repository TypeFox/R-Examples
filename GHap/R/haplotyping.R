#Function: ghap.haplotyping
#License: GPLv3 or later
#Modification date: 5 Feb 2016
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Output haplotype genotype matrix for user-defined haplotype blocks

ghap.haplotyping<-function(
  phase,
  blocks,
  outfile,
  freq=0.05,
  batchsize=500,
  ncores=1,
  verbose=TRUE
){
  
  #Check if phase is a GHap.phase object
  if(class(phase) != "GHap.phase"){
    stop("Argument phase must be a GHap.phase object.")
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
  
  #Identify activated samples
  ids.in <- which(phase$id.in)
  id <- phase$id[ids.in]
  id <- id[1:length(id) %% 2 == 0]
  pop <- phase$pop[ids.in]
  pop <- pop[1:length(pop) %% 2 == 0]
  ids.n <- length(id)
  
  #Output hapsamples file
  write.table(x = cbind(pop,id),file = hapsamples,quote = FALSE,row.names = FALSE,col.names=FALSE)
  
  #Generate batch index
  id1<-seq(1,nrow(blocks),by=batchsize)
  id2<-(id1+batchsize)-1
  id1<-id1[id2<=nrow(blocks)]
  id2<-id2[id2<=nrow(blocks)]
  id1 <- c(id1,id2[length(id2)]+1)
  id2 <- c(id2,nrow(blocks))
  if(id1[length(id1)] > nrow(blocks)){
    id1 <- id1[-length(id1)]; id2 <- id2[-length(id2)]
  }
  
  #Log message
  if(verbose == TRUE){
    cat("Processing ", nrow(blocks), " blocks in:\n", sep="")
    batch <- table((id2-id1)+1)
    for(i in 1:length(batch)){
      cat(batch[i]," batches of ",names(batch[i]),"\n",sep="")
    }
  }
  
  
  block.iter.FUN<-function(i){
    
    outline<-NULL
    
    #Get block info
    block.info <- blocks[i,c("BLOCK","CHR","BP1","BP2")]
    
    #SNPs in the block
    snps <- which(phase$bp >= .subset2(block.info,3) & phase$bp <= .subset2(block.info,4) & phase$marker.in == TRUE)
    phase.A0 <- phase$A0[snps]
    phase.A1 <- phase$A1[snps]
    
    if(length(snps) > 1){
      
      #Subset block
      block.subset <- phase$phase[snps,ids.in]
      haplotypes <- apply(block.subset,MARGIN = 2, paste, collapse="")
      
      #Haplotype library
      lib <- table(haplotypes)/(2*ids.n)
      lib <- lib[which(lib >= freq)]
      
      #Output genotype matrix
      if(length(lib) >= 1){
        for(j in names(lib)){
          phase.geno <- as.integer(haplotypes[1:length(haplotypes) %% 2 == 0] == j)
          phase.geno <- phase.geno + as.integer(haplotypes[1:length(haplotypes) %% 2 == 1] == j)
          j <- unlist(strsplit(j,""))
          phase.allele <- rep(NA,times=length(j))
          phase.allele[j == "0"] <- phase.A0[j == "0"]
          phase.allele[j == "1"] <- phase.A1[j == "1"]
          phase.allele <- paste(phase.allele,collapse="")
          phase.info<-paste(paste(block.info,collapse=" "),phase.allele,sep=" ")
          phase.geno<-paste(phase.geno,collapse=" ")
          outline <- c(outline,phase.info,phase.geno)
        }
      }
      
    }
    return(outline)
  }
  
  
  #Iterate blocks
  nblocks.done <- 0
  for(i in 1:length(id1)){
    
    hapgenotypes.con  <- file(hapgenotypes, open = "a")  
    hapalleles.con  <- file(hapalleles, open = "a")  
    
    #Compute blocks
    mylines<-unlist(mclapply(FUN = block.iter.FUN,id1[i]:id2[i],mc.cores = ncores))
    #Write batch to files
    writeLines(text = mylines[1:length(mylines) %% 2 == 1],con=hapalleles.con)
    writeLines(text = mylines[1:length(mylines) %% 2 == 0],con=hapgenotypes.con)
    #Log message
    if(verbose == TRUE){
      nblocks.done <- nblocks.done + (id2[i]-id1[i]) + 1
      cat(nblocks.done, "blocks written to file\r")
    }
    
    #Close connections
    close(con = hapalleles.con)
    close(con = hapgenotypes.con)
    
  }
  
  
}

