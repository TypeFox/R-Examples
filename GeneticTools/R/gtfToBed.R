gtfToBed <- function(gtf,output="min"){
 genes <- as.character(gtf$V9)
 temp <- strsplit(genes,";")
 genes <- sapply(temp,'[[',4)
 temp <- strsplit(genes," ")
 genes <- sapply(temp,'[[',3)
 
 chroms <- c(1:22,"X","Y")
 bed <- data.frame(Chr="-AAA",Start=-1,Stop=-1,Name="-AAA")
 for(i in 1:length(chroms)){
  temp <- gtf[gtf$V1==chroms[i],c(4,5)]
  chrLabel <- paste("chr",chroms[i],sep="")
  temp <- cbind(rep(chrLabel,nrow(temp)),temp,genes[gtf$V1==chroms[i]])
  colnames(temp) <- c("Chr","Start","Stop","Name")
  bed <- rbind(bed,temp)
 }
 bed <- bed[-1,]
 bed[,4] <- as.character(bed[,4])
 geneNames <- unique(bed[,4])
 #bedOut <- data.frame(ncol=4,nrow=length(geneNames))
 bedOut <- data.frame(Name="-AAA",Chr="-AAA",Start=-1,Stop=-1,stringsAsFactors=FALSE)
 if(output=="min"){
   
   for(i in 1:length(geneNames)){
    temp <- bed[bed[,4]==geneNames[i],c(1,2,3)]
    bedOut[i,1] <- geneNames[i]
    bedOut[i,2] <- temp[1,1]
    bedOut[i,3] <- min(c(temp[,2],temp[,3]))
    bedOut[i,4] <- max(c(temp[,2],temp[,3]))
   }
 } else {
  bedOut <- bed
 }
 bedOut
}
