read.transpose <- function(file, header=TRUE, sep=',',
                           na.strings='---', ...){
##
## 1.  readLines
##
  Txt <- readLines(file)
  Nr <- length(Txt)
  if(Nr<1){
      warning('\nNo data in file ', file)
      attr(Txt, 'headers') <- character(0)
      attr(Txt, 'footers') <- character(0)
      attr(Txt, 'other') <- character(0)
      attr(Txt, 'summary') <- c(headers=0, footers=0, data=0,
                                other=0)
      return(Txt)
  }
##
## 2.  Split into fields
##
  txt <- gsub('\"', '', Txt)
  Split <- strsplit(txt, sep)
  nFields <- sapply(Split, length)
##
## 3.  headers, footers, desired data, and other
##
  datStart <- min(which(nFields>3))
  hdi <- seq(1, length=datStart-1)
  headers <- txt[hdi]
#
  nonHdi <- datStart:Nr
  footStart <- (datStart+min(which(nFields[nonHdi]<4))-1)
  footers <- txt[footStart:Nr]
#
  dati <- datStart:(footStart-1)
#  dat <- Split[dati]
#  dat <- Split[hd.ft[1]:hd.ft[2]]
#  nv <- nrow(dat)
#
  minFields <- min(nFields[dati])
  nv <- (footStart-datStart)
  dat <- vector('list', length=nv)
  for(i in dati){
      combi <- (nFields[i]-minFields+1)
      Si <- Split[[i]]
      di <- paste(Si[2:(combi+1)], collapse=sep)
      dat[[i-datStart+1]] <- c(Si[1], di, Si[(combi+2):nFields[i]])
  }
##
## 4.  Extract column / variable names
##
  if(header){
      h1 <- sapply(dat, '[', 1)
      h2 <- sapply(dat, '[', 2)
      Dat. <- lapply(dat, '[', -(1:2))
      dat. <- do.call(cbind, Dat.)
      colnames(dat.) <- h2
  } else {
      dat. <- do.call(cbind, dat)
      h1 <- rep('', nv)
      h2 <- h1
  }
##
## 5.  Numbers?
##
  datNA <- (dat. %in% na.strings)
  dat0 <- dat.
  dat.[datNA] <- '0'
  out <- as.numeric(dat.)
  if(!any(is.na(out))){
      out[datNA] <- NA
      attributes(out) <- attributes(dat.)
  } else {
      out <- dat.
  }
##
## 6.  Done
##
  attr(out, 'headers') <- headers
  attr(out, 'footers') <- footers
  attr(out, 'summary') <- c(headers=length(headers),
      footers=length(footers), data=length(dat))
  out
}

