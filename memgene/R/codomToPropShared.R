codomToPropShared <- function(alleles, missingData = c(-98,-99), genind=FALSE) {
    
    if (!require(adegenet)) stop("memgene:  Please install the adegenet package to find the proportion of shared alleles", call.=FALSE)
    
    if ((ncol(alleles) %% 2) != 0) { stop("memgene:  Alleles matrix must have paired columns (i.e. diploid genotypes)", call.=FALSE) }

    toGenind <- data.frame(ID=1:nrow(alleles))
  
    for (i in (1:ncol(alleles))[seq(1,ncol(alleles),2)]) {
        toGenind <- cbind(toGenind, paste(alleles[,i], alleles[,i+1], sep="_"))
    }
    names(toGenind) <- c("ID", names(alleles)[seq(1,ncol(alleles),2)])  
    row.names(toGenind) <- toGenind[,"ID"]
    toGenind <- toGenind[,2:ncol(toGenind)]
  
  
    ## Set any loci with missing alleles coded to NA
    toGenind <- apply(toGenind, 2, function(x) as.character(x))
    for (i in 1:ncol(toGenind)) {
        toGenind[grepl(paste(as.character(c(-98,-99)), collapse="|"), toGenind[,i]),i] <- NA
    }
 
    ## Produce genind object using adegenet
    genindObj <- df2genind(toGenind, sep="_", pop=NULL)
    
    ## Find proportion of shared alleles
    if (genind) {
        return(genindObj)
    }
    else {
        return(1-propShared(genindObj))
    }
    
}