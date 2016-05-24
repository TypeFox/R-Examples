getSNPlocations <- function(genotInfo,annot,th){
  colnames(annot) <- c("Chr","Start","End")
  if(!is.null(th)) if(th==0) th <- NULL
  if(!is.null(th)){
    th <- th * 10^6
    genotInfo[,1] <- as.character(genotInfo[,1])
    chrSNPs <- genotInfo[genotInfo[,1]==as.character(annot[1,1]),]
    chrSNPs[,2] <- as.character(chrSNPs[,2])
    chrSNPs[,4] <- as.numeric(as.character(chrSNPs[,4]))  
    lowSNPs <- chrSNPs[chrSNPs[,4]>(annot$Start-th),]
    SNPs <- lowSNPs[lowSNPs[,4]<(annot$End+th),]
    if(is.null(th)) SNPs <- chrSNPs
    output <- list(SNPloc=SNPs,SNPcol=as.numeric(rownames(SNPs)))
  } else {
    output <- list(SNPloc=genotInfo, SNPcol=as.numeric(rownames(genotInfo)))
  }
  output
}