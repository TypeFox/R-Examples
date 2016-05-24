getEnsgInfo <- function(ensg){
  result <- data.frame(ENSG=".", GeneName=".", GeneFull=".", GeneType=".", Summary=".", PubMedHits=".", stringsAsFactors=FALSE)
  for(i in 1:length(ensg)){
    NCBIans <- scan(paste("http://www.ncbi.nlm.nih.gov/gene/?term=",ensg[i],sep=""), what="raw")

  # Get the official name and full name
    offLoc <- which(grepl("Official",NCBIans))
    offName <- NCBIans[offLoc[1]:(offLoc[1]+4)][4]
    offName <- strsplit(strsplit(offName,"class=\"noline\">")[[1]][2],"<span")[[1]][1]
    offFull <- c()
    index <- 0
    detected <- FALSE
    while(!detected){
      offFull[index+1] <- NCBIans[offLoc[2]+3+index]
      detected <- grepl("<span",offFull[index+1])
      index <- index + 1
    }
    offFull[1] <- strsplit(offFull[1],"<dd>")[[1]][2]
    offFull[length(offFull)] <- strsplit(offFull[length(offFull)],"<span")[[1]][1]
    offFull <- paste(offFull,collapse=" ")
  
  # GeneType  
    geneTypeLoc <- which(grepl("type</dt>",NCBIans))
    geneType <- c()
    index <- 0
    detected <- FALSE
    while(!detected){
      geneType[index+1] <- NCBIans[geneTypeLoc[1]+1+index]
      detected <- grepl("</dd>",geneType[index+1])
      index <- index + 1
    }
    geneType[1] <- strsplit(geneType[1],"<dd>")[[1]][2]
    geneType[length(geneType)] <- strsplit(geneType[length(geneType)],"</dd>")[[1]][1]
    geneType <- paste(geneType,collapse=" ")  
    
  # Summary  
    summaryLoc <- which(grepl("Summary</dt>",NCBIans))
    sumText <- c(".")
    index <- 0
    detected <- FALSE
    if(length(summaryLoc)>0){
      while(!detected){
        sumText[index+1] <- NCBIans[summaryLoc[1]+1+index]
        detected <- grepl("</dd>",sumText[index+1])
        index <- index + 1
      }
      sumText[1] <- strsplit(sumText[1],"<dd>")[[1]][2]
      sumText[length(sumText)] <- strsplit(sumText[length(sumText)],"</dd>")[[1]][1]
      sumText <- paste(sumText,collapse=" ")  
    }

  # PubMedHits
    pubMedLoc <- which(grepl("Pubmed\">See",NCBIans))
    pubMedCites <- NCBIans[pubMedLoc+2]
    pubMedCites <- substr(pubMedCites,2,nchar(pubMedCites)-1)
    #pubMedCites <- as.numeric(pubMedCites)
    
    temp <- data.frame(ENSG=ensg[i], GeneName=offName, GeneFull=offFull, GeneType=geneType, Summary=sumText, PubMedHits=pubMedCites, stringsAsFactors=FALSE)
    result <- rbind(result,temp)
  }
  result <- result[-1,]
  result$PubMedHits <- as.numeric(result$PubMedHits)
  rownames(result) <- 1:nrow(result)
  class(result) <- "ensgInfo"
  result
}