# determines the last day an individual was present with respect to a reference date...
# library(EloRating)
# data(adv); data(advpres)
# advpres[1:8, "a"] <- 0
# x <- elo.seq(winner=adv$winner, loser=adv$loser, Date=adv$Date, presence=advpres)
# ID <- "all"
# refdate <- "2010-01-06"

lastdaypresent <- function(x, ID="all", refdate=NULL) {
  if(class(x)=="elo") pm <- x$pmat else stop("so far 'x' must be of class 'elo'...", call.=FALSE)
  if(is.null(refdate)) refdate <- max(x$truedates)
  
  pm <- pm[x$truedates <= refdate, ]
  
  # check for IDs that were never present and remove them from the presence data...
  # might happen if the refdate is early in the sequence and not yet all imigrants have arrived yet...
  if(0 %in% colSums(pm)) { 
      pids <- names(pm)[colSums(pm)==0] 
      pm <- pm[, colSums(pm)>0]
    } else {
      pids <- NA
    } 
  
  res <- x$truedates[apply(pm, 2, function(z)max(which(z==1)))]
  names(res) <- colnames(pm)
  
  # and add ids that were not present yet (if any)
  if(!is.na(pids[1])) { 
    res <- c(res, rep(NA, length(pids))) 
    names(res) <- c(names(pm), pids)
  }
  
  if(ID!="all") {
    if(ID %in% names(res)) res <- as.Date(as.character(res[ID])) else { res <- NA; warning("ID not found...", call.=FALSE)}
  }
  
  return(res)
}

# lastdaypresent(x, "all", refdate="2010-01-02")
# lastdaypresent(x, "xx", refdate="2010-01-02")
# 

