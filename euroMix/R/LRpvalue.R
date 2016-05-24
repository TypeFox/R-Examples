#Wrapper for the functions readData and LRp. Reads data from files and returns a p-value.
LRpvalue <- function(samplefile, victimfile, suspectfile, freqfile, hp, hd, prD, prC ) {
  
  D <- .readMixtureData(freqfile, samplefile, victimfile, suspectfile)
  
  res <- LRp(D$sampleData, D$victimData, D$suspectData, D$db, hp, hd, prD, prC )
  
  return(list(LR=res$LR,pvalue=res$pvalue))
}

#Reads data from files and prepares data matrices for p-value calculations
.readMixtureData <- function(freqfile, samplefile, victimfile, suspectfile){
  
  #Allele frequencies
  alleleFreq <- read.table(freqfile, sep=",", header=TRUE)
  M <- ncol(alleleFreq)-1
  markerNames <- colnames(alleleFreq)[2:(M+1)]
  db <- numeric()
  for(i in 2:(M+1)){
    alleles <- alleleFreq[!is.na(alleleFreq[,i]),c(1,i)]
    Marker <- rep(markerNames[i-1],nrow(alleles))
    dbNew <- cbind(Marker , alleles)
    colnames(dbNew) <- c("Marker","Allele","Frequency")
    db <- rbind(db,dbNew)
  }
  
  #Sample data
  #With AMEL
  #Read first line of file to check how many columns there are
  #testread <- scan(samplefile, what='character', nlines=1, sep=",",quiet=TRUE)
  #sampleDataAll <- read.table(samplefile,sep=",",header=TRUE,na.strings=c("NA","na","NaN","X","Y"),
  #                            colClasses=c('character','character',rep('numeric',length(testread)-2)))
  #ix <- grep('AMEL',sampleDataAll[,2])
  #if(length(ix) > 0) sampleDataAll <- sampleDataAll[-ix,]
  sampleDataAll <- read.table(samplefile,sep=",",header=TRUE)
  sampleData <- sampleDataAll[,3:ncol(sampleDataAll)]
  rownames(sampleData) <- sampleDataAll[,2]
  
  #Victim profile
  #With AMEL
#   victimDataAll <- read.table(victimfile,sep=",",header=TRUE,na.strings=c("NA","na","NaN","X","Y"),
#                               colClasses=c('character','character','numeric','numeric'))
#   ix <- grep('AMEL',victimDataAll[,2])
#   if(length(ix) > 0) victimDataAll <- victimDataAll[-ix,]
  victimDataAll <- read.table(victimfile,sep=",",header=TRUE)
  victimData <- victimDataAll[,3:ncol(victimDataAll)]
  rownames(victimData) <- victimDataAll[,2]  
  
  #Suspect profile(s)
  #With AMEL marker
  #suspectDataAll <- read.table(suspectfile,sep=",",header=TRUE,na.strings=c("NA","na","NaN","X","Y"),
  #                             colClasses=c('character','character','numeric','numeric'))
  #ix <- grep('AMEL',suspectDataAll[,2])
  #if(length(ix) > 0) suspectDataAll <- suspectDataAll[-ix,]
  #Without AMEL
  suspectDataAll <- read.table(suspectfile,sep=",",header=TRUE)
  suspectData <- split(suspectDataAll[,3:4],suspectDataAll[,1])
  for(i in 1:length(suspectData)) rownames(suspectData[[i]]) <- suspectDataAll[1:nrow(suspectData[[1]]),2]
  
  return(list(db=db, sampleData=sampleData, suspectData=suspectData, victimData=victimData))
}
