# creatematrix 14_10_18

creatematrix <- function(eloobject, daterange=NULL, drawmethod="omit", onlyinteracting=FALSE) {
  # set the date range in case it's not specified...
  if(is.null(daterange[1])) { daterange <- c(min(eloobject$truedates), max(eloobject$truedates)) 
  } else { daterange <- as.Date(daterange) }
  
  # get the sequence
  dataseq <- eloobject$logtable
  dataseq$xdate <- eloobject$truedates[1] - 1 + dataseq$Date
  
  #restrict to date range
  dataseq <- dataseq[which(dataseq$xdate >= daterange[1] & dataseq$xdate <= daterange[2]), ]
  
  # create epmty matrix based on presence
  pmat <- eloobject$pmat[which(eloobject$truedates == daterange[1]):which(eloobject$truedates == daterange[2]), ]
  IDS <- sort(colnames(pmat)[which(colSums(pmat)>0)])
  
  mat <- matrix(ncol=length(IDS), nrow=length(IDS), 0)
  colnames(mat) <- rownames(mat) <- IDS
  mat1 <- mat;
  
  # transform factors into characters...
  dataseq$winner <- as.character(dataseq$winner); dataseq$loser <- as.character(dataseq$loser)
  
  # add decided interactions
  xdata <- dataseq[dataseq$draw==FALSE, ]
  xdata <- table(xdata$winner, xdata$loser)
  mat[rownames(xdata), colnames(xdata)] <- xdata
  
  
  # add ties/draws, but separate depending on how they were specified to be treated (if present in the data)
  if(sum(dataseq$draw)>0) {
    
    xdata <- dataseq[dataseq$draw==TRUE, ]
    xdata <- table(xdata$winner, xdata$loser)
    if(drawmethod=="0.5") {
      xdata <- xdata/2
      mat1[rownames(xdata), colnames(xdata)] <- xdata
      mat1 <- mat1 + t(mat1)
      mat <- mat + mat1
    }
    
    #mat1[rownames(xdata), colnames(xdata)] <- xdata
    
    if(drawmethod=="1") { 
      mat1[rownames(xdata), colnames(xdata)] <- xdata
      mat1 <- mat1 + t(mat1)
      mat <- mat + mat1
    }
    
  }
  
  # if "only interacting" was selected: remove those individuals from the matrix that havent interacted...
  if(onlyinteracting) {
    empty <- as.numeric(which(colSums(mat) + rowSums(mat) == 0))
    if(length(empty) > 0) mat <- mat[-empty, -empty]
  }
  
  return(mat)
  
}

