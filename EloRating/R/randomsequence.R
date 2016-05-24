# default number of IDs is 10 and must be smaller than 2601...
# start date: 2000-01-01
# alphabet = T: winner determined by alphabetical order!
# reversals creates x% opposite winner/loser
# ties: should ties/draws be inserted? if yes, x%
# presence: should random presence matrix be generated as well? 
# if NULL:     all IDs present the entire time
# if not NULL: two numbers specify what proportion of IDs was absent what proportion of time
# e.g. c(0.2, 0.1) 20% percent were absent for on average 10% of time



#nID=6; avgIA=15; alphabet = T; reversals=0.1; ties = 5;startdate=as.Date("2000-01-01"); presence=c(0.2, 0.1)
#nID=IDS; avgIA=IA; presence=pres; reversals=reversals; ties=ties; alphabet=T; startdate=as.Date("2000-01-01")

randomsequence <- function(nID=10, avgIA=20, startdate=as.Date("2000-01-01"), alphabet = T, reversals=0.1, ties = NULL, presence=NULL) {
  
  # number of total interactions generated
  # this number may become smaller if presence matrix is implemented?
  totN <- round(nID * avgIA / 2)
  
  # create fictious dates and times
  dates <- seq.Date(from=startdate,by=1,length.out=totN)
  mins <- paste("0",sample(0:59, totN, replace=T), sep="")
  mins[nchar(mins)==3] <- substr(mins[nchar(mins)==3], 2, 3)
  hours <- paste("0",sample(6:17, totN, replace=T), sep="")
  hours[nchar(hours)==3] <- substr(hours[nchar(hours)==3], 2, 3)
  times <- paste(hours, mins, sep=":"); rm(hours, mins)
  
  
  # create 'nID' combinations of IDs from letters
  # 1 letter  = 26
  # 2 letters = 325
  # 3 letters = 2600
  
  if(nID <= 26) {
    IDs <- sort(sample(letters, nID))
  }
  
  if(nID > 26 & nID <= 325) {
    # make sure "in" does not occur as ID because that causes problems...
    com <- combn(26, 2)[, -177]
    samplecom <- com[, sample(1:ncol(com), nID)]
    IDs <- apply(samplecom, 2, function(x) letters[x])
    IDs <- sort(apply(IDs, 2, function(x)paste(x[1],x[2], sep="")))
  }
  
  if(nID > 325 & nID <= 2601) {
    com <- combn(26, 3)
    samplecom <- com[, sample(1:ncol(com), nID)]
    IDs <- apply(samplecom, 2, function(x) letters[x])
    IDs <- sort(apply(IDs, 2, function(x)paste(x[1],x[2],x[3], sep="")))
  }  
  
  
  # create presence matrix
  # if not otherwise specified it is assumed that all IDs were present at all times
  pmat <- matrix(1, ncol=nID, nrow=length(dates)); colnames(pmat) <- IDs
  
  # if specified in call to function, modify the presence matrix
  if(!is.null(presence)) {
    
    # how many and which IDs affected
    ifelse(round(nID*presence[1])!=0, pID <- sample(IDs, round(nID*presence[1])), pID <- sample(IDs, 1))
    
    # modify presence matrix
    # duration of absence:
    dur <- round(presence[2]*length(dates)); if(dur == 0) dur <- 1
    
    i=pID[1]
    for(i in pID) {
      # random reference date
      rdate <- which(dates == sample(dates, 1))
      
      # assign whether "rdate" is exit or enter
      type <- sample(c("exit", "enter"), 1)
      
      # type exit
      if(type=="exit") {
        # exit date
        exdate <- rdate + dur; if(exdate > (length(dates))) exdate <- length(dates)
        pmat[rdate:exdate, i] <- 0
      }
      
      # type enter
      if(type=="enter") {
        # intro date
        indate <- rdate - dur; if(indate < 1) indate <- 1
        pmat[indate:rdate, i] <- 0
      }
    }
    
    
    pmat[2, ] <- pmat[1, ]
    pmat[nrow(pmat)-1, ] <- pmat[nrow(pmat), ]
  }
  
  # make sure that at least 2 individuals were present each day...
  if(0 %in% rowSums(pmat)) for(i in which(rowSums(pmat)==0)) pmat[i, sample(colnames(pmat), 1)] <- 1
  
  if(1 %in% rowSums(pmat)) for(i in which(rowSums(pmat)==1)) pmat[i, sample(colnames(pmat)[pmat[i, ]==0], 1)] <- 1
  
  # create the actual random interactions
  
  IA_table <- cbind(rep("", totN), rep("", totN))
  
  if(alphabet == FALSE) {
    for(i in 1:totN) {
      IA_table[i, ] <- sample(names(which(pmat[i,]==1)), 2)
    }
  }  
  
  if(alphabet == TRUE) {
    for(i in 1:totN) {
      IA_table[i, ] <- sort(sample(names(which(pmat[i,]==1)), 2))
    }
    # create reversals
    if(reversals > 0 & reversals < 1) {
      revsam <- sample(1:totN, round(totN*reversals)); if(length(revsam)==0) revsam <- sample(1:totN, 1)
      IA_table[revsam, ] <- IA_table[revsam, 2:1]
    }
  }  
  
  # create the final sequence data frame
  xdata <- data.frame(Date=dates, Time=times, winner=IA_table[,1], loser=IA_table[,2], Draw=FALSE)
  
  # add proportion of interactions that ended tied/drawn (if specified)
  if(is.null(ties)==FALSE) {
    tiesam <- sample(1:totN, round(totN*ties))
    xdata$Draw[tiesam] <- TRUE; rm(tiesam)
  }
  
  # add Date column to pmat
  pmat <- data.frame(Date=as.Date(dates), pmat)
  
  xdata <- list(seqdat=xdata, pres=pmat)
  class(xdata) <- "randomsequence"
  
  return(xdata)
  
}




#randomsequence(presence=c(0.7, 0.6))

