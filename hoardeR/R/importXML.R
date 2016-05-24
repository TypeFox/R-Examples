
importXML <- function(fa, folder, idTH = 0.8, verbose=TRUE){
  
  seqNames <- names(fa)
  
  # Get the filenames and statistics about how many XML files we have
  fileList <- list.files(folder)
  seqNamesFolder <- unlist(strsplit(fileList,".xml"))
  foundNames <- sum(is.element(seqNames,seqNamesFolder))
  foundNames2 <- sum(is.element(seqNamesFolder,seqNames))
  if(verbose==TRUE){
    cat("Given amount of sequences in fa object:", length(seqNames),"\n")
    cat("XML files in folder:", length(seqNamesFolder),"\n")
    cat("No. of files to import from requested list/given fa object:", foundNames,"(",foundNames/length(seqNamesFolder)*100,"%)\n")  
    cat("FA object names in XML folder:", foundNames2,"(",foundNames2/length(seqNames)*100,"%)\n")  
  }
  # Adjust the importable files to the available ones
  seqNamesFolder <- seqNamesFolder[is.element(seqNamesFolder,seqNames)]
  result <- list()
  # Set which XML files we want to import
    takeThese <- 1:length(seqNamesFolder)
  
  # Go one by one through the XML files and write the results into the list results
  runningIndex <- 1
  for(i in takeThese){
    # Read in the XML file line by line
    res <- readLines(paste(folder,fileList[i],sep=""))
    # Extract the necessary information
    queryLength <- res[grepl("<BlastOutput_query-len>",res)]
    queryLength <- as.numeric(strsplit(strsplit(queryLength,"-len>")[[1]][2],"</Blast")[[1]][1])
    hitOpen <- which(grepl("<Hit>",res)==TRUE)
    hitClose <- which(grepl("</Hit>",res)==TRUE)
    # Set the temporary variables
    underTH <- TRUE
    hitNo <- 1
    hitDev <- NULL
    hitID <- NULL
    hitLen <- NULL
    hitStart <- NULL
    hitEnd <- NULL
    hitChr <- NULL
    # Extract the information as long as we are below the predefined threshold
    while(underTH & hitNo <= length(hitOpen)){
      tempHit <- res[(hitOpen[hitNo]):(hitClose[hitNo])]
      hitDev.temp <- tempHit[grepl("<Hit_def>",tempHit)]
      hitDev.temp <- strsplit(strsplit(hitDev.temp,"<Hit_def>")[[1]][[2]],"</Hit_def>")[[1]][1]
      hitID.temp <- tempHit[grepl("<Hsp_identity>",tempHit)]
      hitID.temp <- as.numeric(strsplit(strsplit(hitID.temp,"<Hsp_identity>")[[1]][[2]],"</Hsp_identity>")[[1]][1])
      hitLen.temp <- tempHit[grepl("<Hsp_align-len>",tempHit)]
      hitLen.temp <- as.numeric(strsplit(strsplit(hitLen.temp,"<Hsp_align-len>")[[1]][[2]],"</Hsp_align-len>")[[1]][1])
      hitStart.temp <- tempHit[grepl("<Hsp_hit-from>",tempHit)]
      hitStart.temp <- as.numeric(strsplit(strsplit(hitStart.temp,"<Hsp_hit-from>")[[1]][[2]],"</Hsp_hit-from>")[[1]][1])
      hitEnd.temp <- tempHit[grepl("<Hsp_hit-to>",tempHit)]
      hitEnd.temp <- as.numeric(strsplit(strsplit(hitEnd.temp,"<Hsp_hit-to>")[[1]][[2]],"</Hsp_hit-to>")[[1]][1])
      
      
      idRatio <- hitID.temp/queryLength
      lenRatio <- hitLen.temp/queryLength
      
      if(idRatio < idTH) underTH <- FALSE
      
      if(underTH){
        hitDev[hitNo] <- hitDev.temp
        hitID[hitNo] <- hitID.temp
        hitLen[hitNo] <- hitLen.temp
        hitChr[hitNo] <- strsplit(strsplit(hitDev.temp,"chromosome ")[[1]][2],",")[[1]][1]
        hitStart[hitNo] <- hitStart.temp
        hitEnd[hitNo] <- hitEnd.temp
      }
      hitNo <- hitNo + 1
    }
    result[[runningIndex]] <- data.frame(Organism = hitDev,
                                         hitID = hitID,
                                         hitLen = hitLen,
                                         hitChr = hitChr,
                                         hitStart = hitStart,
                                         hitEnd = hitEnd, stringsAsFactors=FALSE)
    # Increment the running index
    runningIndex <- runningIndex + 1
    if(verbose) if(runningIndex%%10==0) cat(runningIndex,"XML files processed.\n")
  }
  names(result) <- seqNamesFolder[takeThese]
  class(result) <- "xmlImport"
  result
}