Dates3to1 <- function(data, YMD=c('Year', 'Month', 'Day')){
## 
## 1.  dateCols
##
  dC <- dateCols(data, YMD)
##
## 2.  drop found 
##
  
  dC. <- unlist(dC)
  if(length(dC.)<1){
    Dat <- data 
  } else {
    Dat0 <- data[-dC.]
    nC <- length(dC)
    Dates <- vector(mode='list', nC)
    names(Dates) <- names(dC)
    for(i in seq(length=nC)){
      dCi <- dC[[i]]
      if(length(dCi)>0){
        dati <- data[dCi]
        Dates[[i]] <- Date3to1(dati)
      }
    }
    Dnull <- sapply(Dates, is.null)
    Dat <- cbind(Dat0, as.data.frame(Dates[!Dnull, drop=FALSE]))
  }
##
## 3.  Done 
##
  Dat 
}
