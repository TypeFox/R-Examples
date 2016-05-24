# library(EloRating)
# data(adv)
# el <- elo.seq(winner=adv$winner, adv$loser, Date=adv$Date)
# interactionmatrix <- creatematrix(el)
# 
# runs=10


randomelo <- function(interactionmatrix, runs=2000) {
  # create a sequence from the matrix
  winner <- c()
  loser <- c()
  for(i in 1:nrow(interactionmatrix)) {
    for(j in 1:ncol(interactionmatrix)) {
      if(interactionmatrix[i, j] > 0) {
        winner <- c(winner, rep(rownames(interactionmatrix)[i], interactionmatrix[i, j]))
        loser <- c(loser, rep(colnames(interactionmatrix)[j], interactionmatrix[i, j]))
      }  
    }
  }
  Date <- seq(as.Date("2000-01-01"), as.Date("2000-01-01")+length(winner)-1, by="day")
  # the starting sequence (which will actually not be used, but only randomized versions of it...)
  xdata <- data.frame(Date, winner, loser)
  rm(i,j,winner, loser, Date)
  
  res <- matrix(ncol=length(unique(c(levels(xdata$winner), levels(xdata$loser)))), nrow=runs, 0)
  colnames(res) <- unique(c(levels(xdata$winner), levels(xdata$loser)))
  
  progbar <- txtProgressBar(min = 0, max = runs, style = 3, char=".")
  
  for(i in 1:runs) {
    tempdata <- xdata[sample(1:nrow(xdata)), ]; rownames(tempdata) <- NULL
    tempres <- elo.seq(tempdata$winner, tempdata$loser, tempdata$Date, progressbar=FALSE, runcheck=FALSE)
    res[i, ] <- extract.elo(tempres)[colnames(res)]
    setTxtProgressBar(progbar, i) 
  }
  
  #boxplot(res)
  #apply(res*(-1), 1, rank)
  outp <- list()
  outp[[1]] <- res
  outp[[2]] <- interactionmatrix
  class(outp) <- "randomelo"
  return(outp)
}
