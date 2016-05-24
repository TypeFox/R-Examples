readFinancialCrisisFiles <- function(files, crisisType=7, ...){
##
## 1.  check
##
#  The following generates a NOTE in R CMD check:
#       No visible binding for global variable 'FinancialCRisisFiles'
#  if(missing(files)){
#      data(FinancialCrisisFiles)
#      files <- FinancialCrisisFiles
#  }
#
#  I get the same note without these lines but with
#  files=FinancialCrisisFiles in the function definition
  if(is.character(files))
      files <- financialCrisisFiles(files)
#
  if(!inherits(files, 'financialCrisisFiles'))
      stop('argument files is not of class financialCrisisFiles')
#
  if(length(crisisType)<1)
      stop('length(crisisType)<1:  no data requested')
  if(!is.numeric(crisisType))
      stop('argument crisisType must be integer between 2 and 8')
#
  if(any(abs(crisisType-round(crisisType))>10*.Machine$double.eps))
      stop('argument crisisType must be interger between 2 and 8')
  crisisType <- as.integer(crisisType)
  if(any((crisisType<1) | 8 < crisisType))
      stop('argument crisisType must be integer between 1 and 8')
#  dots
  dots <- list(...)
  dot0 <- dots
  dots$stringsAsFactors <- FALSE
  fileNames <- names(files)
  fi <- file.info(fileNames)
  fileNotFound <- is.na(fi[, 1])
  if(any(fileNotFound))
      stop(fileNames[fileNotFound][1], ' not found')
  nFiles <- length(fileNames)
  outDF <- (length(crisisType)<2)
  Ctries <- unlist(lapply(files, names))
  nCtries <- length(Ctries)
##
## 2.  for each file
##
#  library(gdata)
  iCtry <- 0
  for(i in 1:nFiles){
      iFile <- files[[i]]
      cat('\n', fileNames[i], ': ')
      iNames <- names(iFile)
      nNms <- length(iNames)
      for(j in 1:nNms){
          iCtry <- iCtry+1
          dots$xls <- fileNames[i]
          cat(iNames[j], '')
          dots$sheet <- iNames[j]
          jCtry <- do.call(gdata::read.xls, dots)
          row0 <- (which(jCtry=='1800')-1)
          if(length(row0)!=1)
              stop('data for 1800 not found in sheet ', iNames[j])
          if(outDF){
              if(iCtry<2){
                  N <- nrow(jCtry)
                  nYrs <- (N-row0)
                  outMat <- matrix(NA, nrow=nYrs, ncol=nCtries+1)
                  dimnames(outMat) <- list(NULL, c('year', Ctries))
                  outMat[, 1] <- as.numeric(jCtry[-(1:row0), 1])
                  out <- as.data.frame(outMat)
              }
              out[, iNames[j]] <- as.numeric(jCtry[-(1:row0), 1+crisisType])
          } else {
              if(iCtry<2) {
                  out <- vector('list', length=nCtries)
                  names(out) <- Ctries
              }
              if(length(crisisType)<7){
                  out[[iNames[j]]] <- jCtry[-(1:row0), c(1, 1+crisisType)]
              } else out[[iNames[j]]] <- jCtry
          }
      }
  }
  cat('\n')
##
## done
##
  out
}

