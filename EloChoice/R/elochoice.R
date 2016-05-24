# elochoice 15_10_01

elochoice <- function(winner, loser, kval=100, startvalue=0, runs=1, normprob=FALSE) {
  # , normprob=TRUE # additional argument to allow distinguishing between logistic and normal approach
  # for now, prefer logistic approach
  # normprob <- FALSE

  winner <- as.character(winner); loser <- as.character(loser)

  slfcts <- 0
  if(sum(winner==loser)>0) {
    slfcts <- sum(winner==loser)
    message("data contained ", sum(winner==loser) ," 'self-contests' (identical winner and loser)\nthese cases were excluded from the sequence!")
    ex <- which(winner==loser)
    winner <- winner[-c(ex)]; loser <- loser[-c(ex)]
  }
  if(length(kval) != 1) warning("k has to be of length 1", call. = FALSE)
  if(length(startvalue) != 1) warning("startvalue has to be of length 1", call. = FALSE)


  allids <- sort(unique(c(winner,loser)))
  startvalues <- rep(startvalue[1], length(allids))
  if(normprob) {
    #pmode <- rep(1, length(allids))
    res <- elointnorm(winner, loser, allids, kval[1], startvalues, runs=runs)
  }else{
    #pmode <- rep(2, length(allids))
    res <- eloint(winner, loser, allids, kval[1], startvalues, runs=runs)
  }

  #res <- eloint(winner, loser, allids, kval[1], startvalues, runs=runs, probmode=pmode)

  ratmat <- res[[1]]
  colnames(ratmat) <- allids

  decmat <- res[[4]]
  upsmat <- res[[3]]
  wgtmat <- res[[5]]


  ov <- matrix(nrow = length(allids), ncol = 2)
  rownames(ov) <- allids
  ov[names(table(winner)), 1] <- table(winner)
  ov[names(table(loser)), 2] <- table(loser)
  ov <- cbind(ov, rowSums(ov, na.rm = T))
  ov[is.na(ov)] <- 0
  colnames(ov) <- c("winner", "loser", "total")
  misc <- c(as.character(kval), length(allids), startvalue, runs, length(winner), slfcts)
  names(misc) <- c("kval", "n_allids", "startval", "runs", "totN", "slf")

  res <- list(ratmat=ratmat, decmat=decmat, upsmat=upsmat, wgtmat=wgtmat, misc=misc, ov=ov, ias=cbind(winner, loser))
  class(res) <- "elochoice"
  return(res)
}
