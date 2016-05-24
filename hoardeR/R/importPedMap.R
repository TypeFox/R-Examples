importPedMap <- function(ped, map=NULL, pedSep="\t", pedHeader=FALSE, genoSep=" ", mapSep="\t", mapHeader=FALSE, na.value="0"){
  # If no extension is given, add it to the path
    pedExtension <- substr(ped,nchar(ped)-3,nchar(ped))
    if(pedExtension!=".ped")  ped <- paste(ped,".ped",sep="")
    
  # If no map file is given, we assume that the filename is similar to the ped file, but just without ped extension
    if(is.null(map)) map <- paste(substr(ped,1,nchar(ped)-3),"map",sep="")
    
  # Import the map file
    mapRaw <- read.table(map, sep=mapSep, header=mapHeader, stringsAsFactors=TRUE)
    colnames(mapRaw) <- c("Chr", "SNP", "V3", "Position")
    
  # Import the ped file
    pedRaw <- read.table(ped, sep=pedSep, header=pedHeader, stringsAsFactors=TRUE)
    pedFam <- pedRaw[,1:6]
    colnames(pedFam) <- c("FamilyID",
                          "IndividualID",
                          "PaternalID",
                          "MaternalID",
                          "Sex",
                          "Phenotype")
  
  # Extract the Genotype Matrix
    genoMatrix <- pedRaw[,-(1:6)]
    colnames(genoMatrix) <- mapRaw[,2]
    rownames(genoMatrix) <- pedFam$IndividualID
  
  # Now get the Alleles from the Genotype Matrix
    getAllele <- function(x, genoSep.Internal=genoSep, missing.Internal=na.value){
      genotypes <- levels(x)
      genotypes <- genotypes[genotypes!=paste(missing.Internal,missing.Internal,sep=genoSep.Internal)]
      alleles <- unique(unlist(strsplit(genotypes,genoSep.Internal)))
      alleles <- alleles[alleles!=missing.Internal]
      if(length(alleles)>2) stop ("Sorry. multiallelic marker are currently not supported.")
      if(length(alleles)==1) alleles <- c(alleles,alleles)
      out <- list(allA = alleles[1], allB = alleles[2], genotypes = paste(genotypes, collapse=";"))
      out
    }
  # Store here the alleles   
    alleles <- list()
  # For some reason apply(genoMatrix,2,getAllele) doesn't work, hence for now just a loop:
    for(i in 1:ncol(genoMatrix)){
      alleles[[i]] <- getAllele(genoMatrix[,i])
    }
  
  # Store results in the map data.frame
    mapRaw$alleleA <- sapply(alleles,"[",1)  
    mapRaw$alleleB <- sapply(alleles,"[",2)
    mapRaw$Genotypes <- sapply(alleles,"[",3)
 
    monomorphic <- 0
  # Now transform the genoMatrix to numeric
    for(i in 1:ncol(genoMatrix)){
      alleles <- unique(c(mapRaw[i,]$alleleA,mapRaw[i,]$alleleB))
    # Monomorphic location
      if(length(alleles)>1){
        geno1 <- paste(alleles[1],alleles[1],sep=genoSep)
        geno2 <- paste(alleles[1],alleles[2],sep=genoSep)
        geno3 <- paste(alleles[2],alleles[2],sep=genoSep)
      } else {
        geno1 <- paste(alleles[1],alleles[1],sep=genoSep)
        geno2 <- paste(alleles[1],na.value,sep=genoSep)
        geno3 <- paste(na.value,na.value,sep=genoSep)
        monomorphic <- monomorphic + 1      
      }
      miss <- paste(na.value,na.value,sep=genoSep)
      levels(genoMatrix[,i]) <- list(geno1=geno1, geno2=geno2,geno3=geno3,miss=miss)
      genoMatrix[,i] <- as.numeric(genoMatrix[,i])-1
    }
    
  # Return the results
    out <- list(map=mapRaw, fam=pedFam, geno=genoMatrix, monomorphic=monomorphic)
    class(out) <- "pedMap"
    out
} 
