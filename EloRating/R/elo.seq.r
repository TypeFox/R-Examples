# elo.seq 14_11_19

elo.seq <- function(winner, loser, Date, draw=NULL, presence=NULL, startvalue=1000, k=100, init="average", iterate=0, runcheck=TRUE, progressbar=TRUE) {
  
  if(runcheck) {
    rc <- seqcheck(winner, loser, Date, draw, presence)
    if(sum(rc$checksum[c("IDcheck", "selfinteractions", "startpresence1", "startpresence2", "endpresence1", "endpresence2", "IDmatch", "IA_presencematch", "presenceentries", "datecol", "length", "continouspres")]) > 0) stop("there appear to be some problems with your data, please consider running 'seqcheck()'\notherwise set runcheck=FALSE")
  }
  
  # IDs as character strings and dates as proper dates
  winner <- as.character(winner)
  loser  <- as.character(loser)
  Date   <- as.Date(as.character(Date))
  # check whether there is a column with draw/tie information
  if(is.null(draw)) {
    draw <- rep(FALSE, length(Date))
  }else{
    draw <- as.logical(as.character(draw))
  }
  
  
  # working with dates as integers (first date in sequence is set to 1) (more convenient?)
  ndat <- as.numeric(Date - min(Date) + 1)
  # keep dates for faster referencing (will not be used here but will be part of the output...)
  truedates <- seq(min(Date), max(Date), by="day")
  # all IDs present in the data
  allids <- unique(c(winner, loser))
  
  # matrix with all dates (rows) and all ids (columns) (no rownames or special date-column set since working with numeric dates ('ndat'))
  mat <- matrix(nrow = max(ndat), ncol=length(allids)); colnames(mat) <- allids
  
  ###################################
  #--- formatting presence input ---#
  #------------ START --------------#
  ###################################
  
  # first case, i.e. if nothing is supplied, assume that all individuals were present at all times
  if(is.null(presence)) {
    pmat <- mat; pmat[,] <- 1
  }else{
    # presence is supplied as presence matrix (i.e. actually a data.frame...)
    if(nrow(mat)==nrow(presence)){
      pmat <- presence[,allids]
    }else{
      if(nrow(presence) < nrow(mat)) {
        stop("")#presence has fewer lines than date range of interaction sequence
      }else{
        if(nrow(presence) > nrow(mat)) {
          if("Date" %in% colnames(presence) == FALSE) {
            stop("") #your presence data goes beyond the interaction date range, but lacks a date column (named 'Date')
          }else{
            # if presence data goes beyond date range, cut the presence data accordingly
            mindate <- min(truedates); maxdate <- max(truedates)
            pmat <- presence[which(as.Date(as.character(presence$Date))==mindate) : which(as.Date(as.character(presence$Date))==maxdate), allids]
          }          
        } 
      }    
    }
  }
  
  # make the function stop if there is anything else than 0 and 1 in the matrix
  if(sum(apply(pmat, 2, function(x)(sum(is.na(x))))) > 0) stop("")
  if(sum(apply(pmat, 2, function(x)(sum(x != 0 & x != 1, na.rm=TRUE))))) stop("")
  ###################################
  #--- formatting presence input ---#
  #------------- END ---------------#
  ###################################
  
  ###########################
  #--- create log tables ---#
  #--------- START ---------#
  ###########################
  
  # create the old log table (in a different layout)
  logtable <- data.frame(Date = ndat, winner, loser, Apre=as.numeric(NA), Bpre=as.numeric(NA), Apost=as.numeric(NA), Bpost=as.numeric(NA), draw=draw)
  
  # recent rating table (discard???)
  # recelo <- matrix(nrow=1, ncol=length(allids), NA); colnames(recelo) <- allids
  
  # temporary elo log (extended recelo)
  tempelo <- matrix(ncol=length(allids), nrow=4, NA, dimnames=list(c("recelo", "present", "firstIA", "firstpres"), allids))
  #  # create a matrix with 1=IA and 0=no IA (i.e. did a given ID interact at least once on a given date)
  #  imat <- mat; imat[, ] <- 0; for(i in 1:length(winner)) { imat[ndat[i], c(winner[i], loser[i])] <- 1 }
  # instead: create matrix with number of observed interactions per day
  nmat <- mat; nmat[, ] <- 0
  # for winners
  tabx <- table(Date, winner)
  nmat[as.character(truedates) %in% rownames(tabx) , colnames(tabx)] <- nmat[as.character(truedates) %in% rownames(tabx) , colnames(tabx)] + as.matrix(tabx)
  # for losers  
  tabx <- table(Date, loser)
  nmat[as.character(truedates) %in% rownames(tabx) , colnames(tabx)] <- nmat[as.character(truedates) %in% rownames(tabx) , colnames(tabx)] + as.matrix(tabx)  
  # fill ratings of those individuals that were present at the beginning (i.e. the date of the first interaction/first day of date range) with the startvalue
  startIDs <- colnames(pmat)[pmat[1, ]==1]
  tempelo["recelo", startIDs] <- startvalue
  tempelo["firstpres", startIDs] <- 1
  for(m in startIDs) {tempelo["firstIA", m] <- ndat[min(which(winner==m | loser==m))]}; rm(m)
  # immigrants and their first presence date and first interaction
  imiIDs <- allids[-c(which(allids %in% startIDs))]
  if(length(imiIDs)>0){
    imiIDs <- matrix(ncol=length(imiIDs), nrow=2, 0, dimnames = list(c("firstIA", "firstpres"), imiIDs))
    for(m in colnames(imiIDs)) {imiIDs[1, m] <- ndat[min(which(winner==m | loser==m))]}; rm(m)
    
    # special case if there is only one imigrant in the data set
    
    if(ncol(imiIDs)==1) {
      imiIDs[2,] <- min(which(pmat[,colnames(imiIDs)]==1))
      tempelo[3:4, colnames(imiIDs)] <- imiIDs
    }
    
    if(ncol(imiIDs)>1) {
      imiIDs[2,] <- apply(pmat[,colnames(imiIDs)], 2, function(x){min(which(x==1))})
      # sorting
      imiIDs <- imiIDs[, names(sort(imiIDs[1, ]))]; imiIDs <- imiIDs[, names(sort(imiIDs[2, ]))]
      tempelo[3:4, colnames(imiIDs)] <- imiIDs    
    }
    
    
  }
  tempelo <- tempelo[, names(sort(tempelo[3, ], na.last = TRUE))]
  tempelo <- tempelo[, names(sort(tempelo[4, ], na.last = TRUE))]
  # reorder matrices (to fit same order of names of tempelo...)
  mat <- mat[, colnames(tempelo)]; pmat <- pmat[, colnames(tempelo)]
  # just checking whether first interaction occurred before first presence...
  #if(init=="bottom" | init=="bottom_low") {
  if(min(tempelo[3, ] - tempelo[4, ])<0){
    xx=paste(names(which((tempelo[3, ] - tempelo[4, ])<0)),collapse=", ")
    stop(c("for ID (", xx,") the first interaction occurred before first presence"))
  }
  #}
  
  ###########################
  #--- first interaction ---#
  #-------- START ----------#
  ###########################
  
  # short versions for some objects (can maybe be removed later...)
  W <- winner[1]; L <- loser[1]; D <- ndat[1]
  # the next two lines could be changed to refer to recelo but they might as well not...
  we <- startvalue # most recent rating of winner, in this case = startingvalue
  le <- startvalue # most recent rating of loser, in this case = startingvalue
  
  # calculate new ratings accounting for whether the interaction ended in a draw
  if(draw[1]) {
    newrat <- e.single(we, le, outcome=0, k)
  }else{
    newrat <- e.single(we, le, outcome=1, k)
  }
  
  # insert them in the matrix
  mat[1, c(W, L)] <- newrat
  # fill the logtable as well
  logtable[1, 4:7] <- c(we, le, newrat)
  # fill the recent ratings
  tempelo["recelo", c(W, L)] <- newrat
  
  ###########################
  #--- first interaction ---#
  #--------- END -----------#
  ###########################  
  
  #############################################
  #--- loop through remaining interactions ---#
  #---------------- START --------------------#
  #############################################
  
  # loop for the remaining interactions seperated by init-type:
  log.entry.bottom=cbind(1,"","")
  if(init=="bottom_low"){
    # alle Tiere die am anfang 1000 bekammen bekommen nach der 1. IA wenn sie nicht interagiert haben das minimum
    id.im.min=names(which(tempelo["recelo", ] == startvalue))
    min.rec.elo <- min(tempelo[1,], na.rm=T)
    tempelo["recelo", id.im.min] <- min.rec.elo
    mat[D, id.im.min]<- min.rec.elo
    log.entry.bottom=cbind(1,"",paste(id.im.min,collapse=", "))
  }
  if(init=="bottom_low" | init=="bottom") {
    # progress bar
    if(progressbar) {
      print("loop 1: Elo calculations")
      #flush.console()
      progbar <- txtProgressBar(min = 0, max = length(winner), style = 3, char=".")
    }
    
    for(i in 2:length(ndat)) {
      if(progressbar) setTxtProgressBar(progbar, i) 
      # grab minimum and who has the minimum
      min.rec.elo <- min(tempelo[1,], na.rm=T)
      # short versions for some objects (can maybe be removed later...)
      W <- winner[i]; L <- loser[i]; D <- ndat[i]
      # update tempelo with immigrants indication
      tempelo[2,][colnames(tempelo)%in%names(pmat)[pmat[D, ]==1]]=D
      # if at least one ID is present for the FIRST time: insert a the respective BOTTOM value
      if(length(which(tempelo["firstpres", ] <= D & is.na(tempelo["recelo",]))) > 0) {
        id.f=names(which(tempelo["firstpres", ] <= D & is.na(tempelo["recelo",])))
        tempelo["recelo", id.f] <- min.rec.elo
        mat[D, id.f]<- min.rec.elo
      }
      if(!exists("id.f")){id.f=""}
      # current ratings
      we <- tempelo["recelo", W]; le <- tempelo["recelo", L]
      # calculate new ratings
      #newrat <- e.single(we, le, outcome=1, k)
      
      # calculate new ratings accounting for whether the interaction ended in a draw
      if(draw[i]) {
        newrat <- e.single(we, le, outcome=0, k)
      }else{
        newrat <- e.single(we, le, outcome=1, k)
      }
      
      # insert them in the matrix
      mat[D, c(W, L)] <- newrat
      # insert them in the recentelo
      tempelo["recelo", c(W, L)] <- newrat
      # fill the logtable as well
      logtable[i, 4:7] <- c(we, le, newrat)
      # all noninteracting animals and previously lowest rankers get the lowest updated lowest rank
      who.min.rec.elo <- names(which(tempelo[1,]==min.rec.elo))
      min.rec.elo <- min(tempelo[1,], na.rm=T)
      id.min=setdiff(who.min.rec.elo,c(W, L))
      tempelo[1,id.min]<- min.rec.elo
      mat[D, id.min]<- min.rec.elo
      if(!exists("id.min")){id.min=""}
      log.entry.bottom=c(log.entry.bottom, cbind(D,paste(id.f,collapse=", "),paste(id.min,collapse=", ")))
      rm(id.f,id.min)
    }
    log.entry.bottom=as.data.frame(matrix(log.entry.bottom, ncol=3,nrow=(length(ndat)), byrow = T))
    names(log.entry.bottom)=c("Date", "new.entry", "set.to.bottom.value")
  }
  # 'bottom section' END
  
  if(init=="average") {
    # progress bar
    if(progressbar) {
      print("loop 1: Elo calculations")
      #flush.console()
      progbar <- txtProgressBar(min = 0, max = length(winner), style = 3, char=".")
    }
    
    log.entry.bottom=cbind(1,"")
    for(i in 2:length(ndat)) {
      if(progressbar) setTxtProgressBar(progbar, i) 
      # short versions for some objects (can maybe be removed later...)
      W <- winner[i]; L <- loser[i]; D <- ndat[i]
      # update tempelo with immigrants indication
      tempelo[2,][colnames(tempelo)%in%names(pmat)[pmat[D, ]==1]]=D
      # if at least one ID is present for the FIRST time: insert a the respective AVERAGE value
      if(length(which(tempelo["firstpres", ] <= D & is.na(tempelo["recelo",]))) > 0) {
        avg.elo <- round(mean(tempelo["recelo", tempelo["present", ]<= D], na.rm=T))
        id.f=names(which(tempelo["firstpres", ] <= D & is.na(tempelo["recelo",])))
        tempelo["recelo", id.f] <- avg.elo
        mat[D, id.f]<- avg.elo
      }
      if(!exists("id.f")){id.f=""}
      # current ratings
      we <- tempelo["recelo", W]; le <- tempelo["recelo", L]    
      # calculate new ratings
      #newrat <- e.single(we, le, outcome=1, k)
      
      # calculate new ratings accounting for whether the interaction ended in a draw
      if(draw[i]) {
        newrat <- e.single(we, le, outcome=0, k)
      }else{
        newrat <- e.single(we, le, outcome=1, k)
      }
      
      # insert them in the matrix
      mat[D, c(W, L)] <- newrat
      # insert them in the recentelo
      tempelo["recelo", c(W, L)] <- newrat  
      # fill the logtable as well
      logtable[i, 4:7] <- c(we, le, newrat)
      log.entry.bottom=c(log.entry.bottom, cbind(D,paste(id.f,collapse=", ")))
      rm(id.f)
    }
    log.entry.bottom=as.data.frame(matrix(log.entry.bottom, ncol=2,nrow=(length(ndat)), byrow = T))
    names(log.entry.bottom)=c("Date", "new.entry")
  }
  # 'average section' END
  
  if(progressbar) close(progbar)
  #############################################
  #--- loop through remaining interactions ---#
  #----------------- END ---------------------#
  #############################################
  
  ################################
  #--- matrix fill LARS style ---#
  #---------- START -------------#
  ################################
  
  # start with the original rating matrix
  lmat <- rbind(rep(NA,ncol(mat)),mat)
  lmat[1, startIDs] <- startvalue
  # fill ratings with an ID's rating from the day before if on a given day it is NA
  for(i in 2:nrow(lmat)){
    lmat[i, is.na(lmat[i, ])] <- lmat[i-1, is.na(lmat[i, ])]
  }
  # remove ratings on days which IDs were not present (based on presence data)
  lmat=lmat[-1,]
  lmat[pmat==0] <- NA
  
  ################################
  #--- matrix fill LARS style ---#
  #----------- END --------------#
  ################################
  
  ####################################
  #--- matrix fill CHRISTOF style ---#
  #------------ START ---------------#
  ####################################
  
  # start with the original rating matrix
  cmat <- mat; 
  # first: for those IDs that were present at the beginning (startIDs): backfill the first rating back up to day 1
  needtobefilled <- names(which(is.na(mat[1, startIDs])))
  if (length(needtobefilled)>0) {
    for(ID in needtobefilled) { cmat[1, ID] <- lmat[min(which(complete.cases(lmat[, ID]))), ID]}
  }
  
  #  # fill FIRST rating (the day of entry) if no observation was observed that day (similar to above for lmat)
  #  # only first is filled, as interpolation follows
  #  for(C in 1:ncol(mat)) {
  #    cbind(cmat[,C], pmat[,C])
  #    (firstelo <- min(which(complete.cases(cmat[, C]))))
  #    (firstpresent <- min(which(pmat[, C]==1)))
  #    if(firstpresent < firstelo) {
  #      if(init=="average") cmat[firstpresent, C] <- round(mean(lmat[firstpresent, names(which(pmat[firstpresent, ]==1))], na.rm=TRUE))
  #      if(init=="bottom")  cmat[firstpresent, C] <-        min(lmat[firstpresent, names(which(pmat[firstpresent, ]==1))], na.rm=TRUE)
  #    }
  #    cmat[, C]
  #  }
  
  # which IDs were observed only at first day
  only.f <- as.numeric(which(apply(mat, 2, function(x) (length(unique(x))))==2 & apply(mat, 2, function(x) (!is.na(x[1])))))
  # fill ratings with an ID's rating from the day before if on a given day it is NA
  if(length(only.f)>0){
    for(i in only.f){
      for(j in 2:nrow(cmat)){
        cmat[j, i] <- cmat[j-1, i]
      }
    }
  }
  # which IDs were observed (or 'guessed') at least twice
  morethan1 <- as.numeric(which(apply(lmat, 2, function(x) (length(unique(x)))) > 2) )
  for(i in morethan1) { cmat[,i] <- round(na.approx(cmat[,i], rule = 2)) }
  #cmat[, morethan1] <- apply(cmat[, morethan1], 2, function(x) round(na.approx(x, rule = 2)))
  cmat[pmat==0] <- NA
  
  ####################################
  #--- matrix fill CHRISTOF style ---#
  #------------- END ----------------#
  ####################################
  
  ################################
  #--- stability calculations ---#
  #----------- START ------------#
  ################################
  
  # subfunction to calculate standardized ratings, used for weighing the stability index
  # standardized ratings range between 0 (for the individual with the lowest rating) and 1 (for the individual with the highest rating on the given day)
  elo.stdz <- function(ratings) {
    ratings <- ratings - min(ratings, na.rm=TRUE)
    return(ratings/max(ratings, na.rm=TRUE))
  }
  
  # get the IDs of all individuals present in the data (and the total number of IDs)
  # AllIDs <- colnames(drm); nIDs <- length(AllIDs)
  # create empty vectors for the three variables of interest
  rankdiffs <- c(); Idspresent <- c(); eloweights <- c()
  
  # progress bar
  if(progressbar) {
    print("loop 2: Stability calculations")
    #flush.console()
    progbar <- txtProgressBar(min = 0, max = nrow(cmat), style = 3, char=".")
  }
  
  # this loop calculates Ci for each day (except for the first one)
  for(u in 2:nrow(cmat)) {
    if(progressbar) setTxtProgressBar(progbar, u) 
    # calculates the ranks the day before the actual day
    r1 <- rank(cmat[u-1, ] * (-1), na.last = NA, ties.method = c("average"))
    # calculates the ranks on the test day
    r2 <- rank(cmat[u, ]   * (-1), na.last = NA, ties.method = c("average"))
    # which IDs were present on both days
    present <- c(names(r1), names(r2))[duplicated(c(names(r1), names(r2)))]  
    # if one animal leaves, the index increases the ranks of all individuals below, i.e. if no other rank change occurs, the rankdifference will be zero in such a case
    # if(length(r1) > length(r2)) { # faulty line!!!!
    if(length(which(!names(r1) %in% names(r2))) > 0) {
      #leavers <- names(which(table(c(names(r1), names(r2))) == 1)) # faulty as well
      leavers <- names(r1)[which(!names(r1) %in% names(r2))]
      for(n in 1:length(leavers)) {
        r1[which(r1 > r1[leavers[n]])] <- r1[which(r1 > r1[leavers[n]])] - 1
      }
      r1 <- r1[-c(which(names(r1) %in% leavers))]
      rm(leavers)
    }  
    # calculate the weights of change (if there is none, the weight is '0')
    standardratings <- elo.stdz(cmat[u - 1, present])
    changers <- r1[r1[present] != r2[present]]
    stabweight <- 0
    if(length(changers) > 0) {
      stabweight <- as.numeric(standardratings[names(changers)[changers==min(changers)][1]])
      rm(changers)
    }
    
    # calculate the sum of the absolute differences in the two rankings
    rankdiffs <- c(rankdiffs, sum(abs(r2[present] - r1[present])))
    # how many individuals were present on both days
    Idspresent <- c(Idspresent, length(present))
    # the standardized elo rating of the highest rated individual involved in a rank change
    eloweights <- c(eloweights, stabweight)
    rm(stabweight, present, r2, r1)  
  } # end of loop through dailyratingmatrix ('Christof - style')
  
  dte <- seq(from=as.Date(min(Date)), to=as.Date(max(Date)), by=1)
  stability <- data.frame(date=dte[2:length(dte)], Idspresent, rankdiffs, eloweights)
  
  if(progressbar) close(progbar)
  
  ################################
  #--- stability calculations ---#
  #------------ END -------------#
  ################################
  
  ###########################
  #--- final log data... ---#
  #-------- START ----------#
  ###########################
  
  # get some more 'log' data...
  
  
  
  misc <- matrix(c("init",init,
                   "k", k, 
                   "startvalue", startvalue,
                   "minDate", as.character(as.Date(min(Date))), 
                   "maxDate", as.character(as.Date(max(Date))),
                   "nID", length(table(c(winner,loser))),
                   "IAmin", min(table(c(winner,loser))),
                   "IAmax", max(table(c(winner,loser))),
                   "IAmean", round(mean(table(c(winner,loser))),1),
                   "IAmedian", round(median(table(c(winner,loser))),1),
                   "nIA", length(winner),
                   "IAperDay", round(mean(rowSums(nmat)/2),1),
                   "draws", round(sum(draw)/length(draw),2)
  ),ncol=2,byrow=T)
  rownames(misc) <- misc[,1]; misc <- misc[,2]
  
  ###########################
  #--- final log data... ---#
  #--------- END -----------#
  ###########################
  
  if(length(imiIDs)==0){
    logtable=logtable
  }else{
    if(init=="average"){
      logtable=data.frame(logtable,new.entry=log.entry.bottom[,-1])
    }else{
      logtable=data.frame(logtable, log.entry.bottom[,-1])
    }
  }
  # return and 'class' results
  res <- list(mat=mat, lmat=lmat, cmat=cmat, pmat=pmat, nmat=nmat, logtable=logtable, stability=stability, truedates=truedates, misc=misc, allids=sort(allids))
  
  class(res) <- "elo"
  return(res)
}

