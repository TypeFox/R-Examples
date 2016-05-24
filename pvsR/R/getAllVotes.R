##' Get several votes
##' 
##' This function is essentially a  wrapper around Votes.getBillActionVotes() specified for large amount of requests.
##' @usage getAllVotes(actionId, batchsize=100, pause=0, backupfile="votes.list.Rdata")
##' @param actionId a character string or list of character strings with the action ID(s) (see references for details)
##' @param batchsize numerical, indicating how many actionIds should be processed in one batch (defaults to 100).
##' @param pause numerical, indicating how long (in seconds) the download process should be paused after each batch (defaults to 0)
##' @param backupfile character string for the path/file-name of the Rdata-file where the data should be saved (batch-wise) during the download process (default: "votes.list.Rdata").
##' @return A data frame with a row for each vote and columns with the following variables describing the vote:\cr votes.vote*.candidateId,\cr votes.vote*.candidateName,\cr votes.vote*.officeParties,\cr votes.vote*.action.
##' @details This functions splits large requests into several batches. The requests are then processed batch-wise and are saved on the local disc to make sure that not too much RAM is assigned to the pvsR task.
##' @references http://api.votesmart.org/docs/Votes.html\cr
##' Use Votes.getBill() or Votes.getByOfficial() to get a list of action IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get all officials of a certain state
##' \dontrun{bill <- Votes.getBill("17623", separate=c("actions", "sponsors"))}
##' \dontrun{actionids <- bill$actions$actionId}
##' # get all votes on this acti
##' \dontrun{votes <- getAllVotes(actionids, batchsize=2)}
##' \dontrun{head(votes)}
##' @export


getAllVotes <-
  function(actionId, batchsize=100, pause=0, backupfile="votes.list.Rdata") {
    
    n <- length(actionId)
    rest <- n%%batchsize
    
    chunks.upper <- seq(from = batchsize, to = n, by = batchsize)
    
    
    if (rest != 0) {
      
      chunks.upper[length(chunks.upper) + 1] <- chunks.upper[length(chunks.upper)] + rest
      
    }
    
    chunks.lower <- c(1,chunks.upper[-length(chunks.upper)] + 1)
    
    
    
    # prepare for loop over all chunks
    chunks <- data.frame(lower=chunks.lower, upper=chunks.upper)
    pb <- txtProgressBar(min = 0, max = nrow(chunks), style = 3)
    
    votes.list <- as.list(1:nrow(chunks))
    save(votes.list, file=backupfile) # to be saved and loaded in each loop
    
    # process queries chunkwise
    for (i in 1:nrow(chunks)) {
      
      Sys.sleep(pause)
      
      first <- chunks$lower[i]
      last <- chunks$upper[i]
      
      
      cIds <- actionId[first:last]
      
      votes <- Votes.getBillActionVotes(cIds)
      
      
      load(backupfile)
      votes.list[[i]] <- votes
      save(votes.list, file=backupfile)
      rm(votes.list )
      gc(, verbose=FALSE) # clean memory
      
      setTxtProgressBar(pb, i)
      
    }
    
    load(backupfile)
    allvotes <- dfList(votes.list)
    
    return(allvotes)
    
    
  }