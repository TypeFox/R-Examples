# seqcheck 14_11_19


seqcheck <- function(winner, loser, Date, draw=NULL, presence=NULL) {
  
  # Part 1:
  # check some general issues 
  
  # creating checksum
  checksum <- rep(0, 14)
  names(checksum) <- c("IDcheck", "selfinteractions", "presence", "startpresence1", "startpresence2", "endpresence1", "endpresence2", "IDmatch", "IA_presencematch", "presenceentries", "datecol", "length", "singledayobs", "continouspres")
  
  Date <- as.Date(as.character(Date))
  
  # check whether all vectors contain same number of entries
  if(length(winner)!=length(loser) | length(winner)!=length(Date)) checksum["length"] <- 1
  
  # check draw/tie vector
  if(is.null(draw)) draw <- rep(FALSE, length(winner))
  if(length(winner)!=length(draw)) checksum["length"] <- 1
  
  
  # the remaining checks are conditional on the fact that length of vectors match...
  
  if(checksum["length"]==0) {
    
    datasequence <- data.frame(Winner=winner, Loser=loser, Date=Date)
    
    # check for integrity of IDs 
    winners <- as.character(datasequence[, "Winner"])
    losers <- as.character(datasequence[, "Loser"])
    allIDs <- sort(unique(c(winners, losers)))
    if(length(allIDs)==length(unique(tolower(allIDs)))) {
      IDcheck <- "everything seems alright with capitalization of IDs"
    }else{
      IDcheck <- "seems to be a problem with capitalization of IDs?"
      checksum[1] <- 1
    }
    
    
    # check whether IDs had interactions with themselves...
    selfinteractions <- paste("There is (are) ", length(which(winners==losers)), " case(s) in which loser ID equals winner ID", sep="")
    if(length(which(winners==losers)) > 0) checksum[2] <- 1
    
    # check whether presence data is given
    if(!is.null(presence)) {
      presenceD <- "presence data is supplied"
    }else{
      presenceD <- "presence data is not supplied"
      checksum["presence"] <- 1
    }
    
    # check whether there is a column named Date in the presence matrix
    if(!is.null(presence)) {
      if("Date" %in% colnames(presence)) {
        datecol <- "Date column found"
      }else{
        datecol <- "no 'Date' column found in supplied presence data"
        checksum["datecol"] <- 1
      }
    }
    
    # check whether there are gaps in the presence data...
    if(!is.null(presence) & checksum["datecol"]==0) {
      if(nrow(presence) < as.numeric(diff(range(presence$Date)))+1) {
        checksum["continouspres"] <- 1
        continouspres <- "there appear to be gaps in your presence data (missing days?)"
      }
    } else {
      continouspres <- "not checked"
    }
    
    
    
    # check whether date range in presence is the same as in sequence data 
    START <- NA
    END <- NA
    if(!is.null(presence) & checksum["datecol"]==0) {
      DATESpres <- as.Date(as.character(presence[,which(colnames(presence) %in% c("Date", "date"))]))
      DATESdata <- unique(as.Date(as.character(datasequence[,which(colnames(datasequence) %in% c("Date", "date"))])))
      if(min(DATESpres)<min(DATESdata))  {START <- "presence starts earlier than data"; checksum[4] <- 1} # actually, not a problem per se
      if(min(DATESpres)>min(DATESdata))  {START <- "presence starts AFTER data -> PROBLEM!"; checksum[5] <- 1}
      if(min(DATESpres)==min(DATESdata)) {START <- "presence starts at the same date than data -> GOOD!"}
      if(max(DATESpres)<max(DATESdata))  {END   <- "presence stops BEFORE data -> PROBLEM!"; checksum[6] <- 1}
      if(max(DATESpres)>max(DATESdata))  {END   <- "presence continues beyond data"; checksum[7] <- 1} # actually, not a problem per se
      if(max(DATESpres)==max(DATESdata)) {END  <- "presence stops at the same date than data -> GOOD!"}
    }
    
    # check whether IDs match in data and presence 
    if(!is.null(presence)) {  
      IDdata <- sort(allIDs)
      IDpres <- sort(names(presence[,2:ncol(presence)]))
      
      IDmatch <- "IDs in presence and data match -> GOOD!"
      
      # check whether 
      if(length(which(!IDpres %in% IDdata)) > 0) { 
        IDmatch1 <- IDpres[which(!IDpres %in% IDdata)] 
      }else{
        IDmatch1 <- "none"
      }
      
      if(length(which(!IDdata %in% IDpres)) > 0) { 
        IDmatch2 <- IDdata[which(!IDdata %in% IDpres)] 
      }else{
        IDmatch2 <- "none"
      }
      
      if(IDmatch1[1]!="none" | IDmatch2[1]!="none") {
        IDmatch <- c(paste("the following IDs occur in the presence data but NOT in the data sequence:", IDmatch1), paste("the following IDs occur in the data sequence but NOT in the presence data:", IDmatch2))
        checksum[8] <- 1
      }
      
    }
    
    # check whether IDs were actually present on the dates of their interactions
    # note that IDs that occur in the data sequence BUT NOT in the presence are ignored here!
    IA_presencematch  <- NA
    IA_presencematchN <- NA
    
    nmatrows <- length(seq(from=min(as.Date(Date)), to=max(as.Date(Date)), by="day"))
    if(!is.null(presence) & checksum["datecol"]==1) checksum[9] <- 2
    
    if(!is.null(presence) & checksum["datecol"]==0 ) {
      if(nmatrows==nrow(presence)) {
        IA_presencematch <- c()
        for(i in 1:nrow(datasequence)) {
          if(sum(!(IDmatch2 %in% c(winners[i],losers[i]))) > 0 & IDmatch2[1]=="none") {
            if(sum(presence[which(datasequence[i, "Date"] == presence[, "Date"]) , c(winners[i],losers[i])], na.rm=TRUE) != 2) {
              IA_presencematch <- c(IA_presencematch, i)
            }
          }
        }
        if(is.null(IA_presencematch)) {
          IA_presencematch <- "all IDs were present on their interaction dates"
          IA_presencematchN <- 0
        }else{
          IA_presencematchN <- IA_presencematch
          IA_presencematch <- paste("during", length(IA_presencematchN), "interactions, IDs were absent according to the presence data:", sep=" ")
          checksum[9] <- 1
        }
      }
    }
    
    # check whether there are entries other than 0's and 1's in the presence
    if(!is.null(presence)) {
      temp <- as.numeric(apply(presence[,2:ncol(presence)], 2, function(xx)length(which(xx==0 | xx==1))))
      presenceentries <- "all presence entries are either 1 or 0 --> GOOD"
      if(length(which(temp!=nrow(presence))) > 0)  {
        presenceentries <- "at least one presence entry is not 1 or 0 --> PROBLEM"
        checksum[10] <- 1                                        
      }
    }  
    
    # check for cases in which IDs were observed only on a single day (even though multiple times is possible, but that doesnt make a difference...)
    temp <- rbind(rowSums(table(winner, Date)>0)[allIDs], rowSums(table(loser, Date)>0)[allIDs])
    colnames(temp) <- allIDs
    sIDs <- "none"
    if(1 %in% colSums(temp, na.rm=T)) {
      checksum["singledayobs"] <- 1
      sIDs <- colnames(temp)[colSums(temp, na.rm=T)==1]
    }
    
    
    if(!is.null(presence)) {
      res <- list(checksum = checksum, IDcheck=IDcheck, selfinteractions=selfinteractions, presence=presenceD, startpresence = START, endpresence = END, IDmatch=IDmatch, IA_presencematch = IA_presencematch, IA_presencematchN = IA_presencematchN, presenceentries = presenceentries, IDmatch1 = IDmatch1, IDmatch2 = IDmatch2, datecol = datecol, singledaycases=sIDs)
      class(res) <- "sequencecheck"
    }
    if(is.null(presence)) { 
      res <- list(checksum = checksum, IDcheck=IDcheck, selfinteractions=selfinteractions, presence=presenceD, singledaycases=sIDs, continouspres=continouspres)
      class(res) <- "seqchecknopres"
    }
    
    
  } # end part (conditional on vector length match)
  
  if(checksum["length"]==1) { 
    res <- list(checksum = checksum, IDcheck=NA, selfinteractions=NA, presence=NA)
    
    class(res) <- "seqchecknopres"
    
  }
  
  return(res)
}




