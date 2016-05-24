##' Fetch data on all local (city- or county-) officials
##' 
##' This function is essentially a  wrapper around Local.getOfficials().
##' @usage getAllLocalOfficials(locality="counties", batchsize=50,
##'  pause=0, backupfile="locofs.list.Rdata")
##' @param locality a character string indicating whether data on county-officials ("counties") or city-officials ("cities") should be downloaded.
##' @param batchsize numerical, indicating the number of requests that should be processed in one batch (defaults to 50).
##' @param pause numerical, indicating how long (in seconds) the download process should be paused after each batch (defaults to 0)
##' @param backupfile character string for the path/file-name of the Rdata-file where the data should be saved (batch-wise) during the download process (default: "locofs.list.Rdata").
##' @return A data frame with a row for each official and columns with the following variables describing the official:\cr candidatelist.candidate*.candidateId,\cr candidatelist.candidate*.firstName,\cr candidatelist.candidate*.nickName,\cr candidatelist.candidate*.middleName,\cr candidatelist.candidate*.lastName,\cr candidatelist.candidate*.suffix,\cr candidatelist.candidate*.title,\cr candidatelist.candidate*.electionParties,\cr candidatelist.candidate*.electionDistrictId,\cr candidatelist.candidate*.electionStateId,\cr candidatelist.candidate*.officeParties,\cr candidatelist.candidate*.officeDistrictId,\cr candidatelist.candidate*.officeDistrictName,\cr candidatelist.candidate*.officeStateId,\cr candidatelist.candidate*.officeId,\cr candidatelist.candidate*.officeName,\cr candidatelist.candidate*.officeTypeId.
##' @details This functions splits large requests into several batches. The requests are then processed batch-wise and are saved on the local disc to make sure that not too much RAM is assigned to the pvsR task.
##' @references http://api.votesmart.org/docs/Local.html\cr
##' Use Local.getCounties() or Local.getCities() to get a list of local IDs. 
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of all local officials
##' \dontrun{all_countyofficials <- getAllLocalOfficials()}
##' \dontrun{head(officials)}
##' @export


getAllLocalOfficials <-
  function(locality="counties", batchsize=50, pause=0, backupfile="locofs.list.Rdata") {
    
    if(locality=="counties"){
      
      counties <- getAllCounties()
      localId <- counties$localId
      
    }
    
    if(locality=="cities"){
      
      counties <- getAllCities()
      localId <- counties$localId
      
    }
    
    n <- length(localId)
    rest <- n%%batchsize
    
    chunks.upper <- seq(from = batchsize, to = n, by = batchsize)
    
    
    if (rest != 0) {
      
      chunks.upper[length(chunks.upper) + 1] <- chunks.upper[length(chunks.upper)] + rest
      
    }
    
    chunks.lower <- c(1,chunks.upper[-length(chunks.upper)] + 1)
    
    
    
    # prepare for loop over all chunks
    chunks <- data.frame(lower=chunks.lower, upper=chunks.upper)
    pb <- txtProgressBar(min = 0, max = nrow(chunks), style = 3)
    
    ofs.list <- as.list(1:nrow(chunks))
    save(ofs.list, file=backupfile) # to be saved and loaded in each loop
    
    # process queries chunkwise
    for (i in 1:nrow(chunks)) {
      
      Sys.sleep(pause)
      
      first <- chunks$lower[i]
      last <- chunks$upper[i]
      
      
      locIds <- localId[first:last]
      
      bios <- Local.getOfficials(locIds)
      
      
      load(backupfile)
      ofs.list[[i]] <- bios
      save(ofs.list, file=backupfile)
      rm(ofs.list )
      gc(, verbose=FALSE) # clean memory
      
      setTxtProgressBar(pb, i)
      
    }
    
    load(backupfile)
    allLocalOffs <- dfList(ofs.list)
    
    return(allLocalOffs)
    
  }
