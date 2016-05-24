readNIPA <- function(files, sep.footnote='/', ...){
##
## 1.  read.transpose
##
  nFiles <- length(files)
  Data <- vector('list', length=nFiles)
  for(i in 1:nFiles){
      Data[[i]] <- read.transpose(files[i], ...)
  }
#  Data <- lapply(files, read.transpose)
  isNum <- sapply(Data, is.numeric)
  if(!all(isNum)){
      NonNum <- sum(!isNum)
      badFiles <- paste(files[!isNum], collapse=', ')
      oops <- paste('Nonnumerics found in file', c('', 's')[1+(NonNum>1)],
                    badFiles, sep='')
      stop(oops)
  }
##
## 2.  Summary stats
##
  sumStats <- lapply(Data, attr, 'summary')
  SumStats <- do.call(rbind, sumStats)
  files. <- strsplit(files, '/')
  fileNms <- sapply(files., function(x)x[length(x)])
  fileNm. <- sub('.csv', '', fileNms)
  rownames(SumStats) <- fileNm.
  year1 <- sapply(Data, function(x)min(x[, 1]))
  yearn <- sapply(Data, function(x)max(x[, 1]))
  overlap <- c(NA, yearn[-nFiles]-year1[-1]+1)
#
  lastRev <- function(x){
      hd <- attr(x, 'header')
      upd <- grep('Last Rev', hd, value=TRUE)
      Upd <- sub('Last Revised on:', '', upd)
      UPD <- strsplit(Upd, '-')[[1]][1]

      lastRev <- as.Date(UPD, format=' %B %d,%Y')
  }
  lastUpdate <- lapply(Data, lastRev)
  lastRev. <- lastUpdate[[1]]
  for(i in seq(2, length=nFiles-1))
      lastRev.[i] <- lastUpdate[[i]]
  sumSt <- cbind(as.data.frame(SumStats),
                 yearFirst=year1, yearLast=yearn, yearsOverlap=overlap,
                 lastRevision=lastRev.)
##
## 3.  common initial colnames
##
  maxVars <- min(SumStats[, 3])
  varNms <- lapply(Data, colnames)
  VarNms <- lapply(varNms, function(x){
      vn <- strsplit(x, sep.footnote)
      sapply(vn, '[', 1)
  } )
  varNms1 <- VarNms[[1]]
  for(i in seq(2, length=nFiles-1)){
      iVars <- 1:maxVars
      comp <- (varNms1[iVars]==VarNms[[i]][iVars])
      if(comp[1]){
          maxVars <- max(which(comp))
      } else stop('no common variables between the first and file ', i)
  }
  keepVars <- 1:maxVars
  keepNames <- varNms1[keepVars]
##
## 4.  rbind
##
  rmsd <- rep(NA, nFiles)
  if(nFiles>1){
#      nrow1 <- with(sumSt, yearLast[1]-yearFirst[1]-overlap[2]+1)
      Dat1 <- Data[[1]][, keepVars]
#      Dat1 <- Dati[1:nrow1, ]
#
      for(i in 2:nFiles){
#          Datm1 <- Dat1
          Dati <- Data[[i]][, keepVars]
          if(overlap[i]>0){
              n1 <- nrow(Dat1)
              n0 <- (n1-overlap[i])
              O1 <- Dat1[n0+1:overlap[i], ]
              Oi <- Dati[1:overlap[i], ]
              na1 <- is.na(O1)
              nai <- is.na(Oi)
              na.1.i <- !(na1 | nai)
              D1 <- ((Oi[na.1.i]-O1[na.1.i])/O1[na.1.i])
              rmsd[i] <- sqrt(mean(D1^2))
              O1[na.1.i] <- (O1[na.1.i]+Oi[na.1.i])/2
              na1.i <- (na1 & !nai)
              O1[na1.i] <- Oi[na1.i]
              Dat1 <- rbind(Dat1[1:n0, ], O1)
          }
          Dat1 <- rbind(Dat1, Dati[(overlap[i]+1):nrow(Dati), ])
      }
  } else Dat1 <- Data[[1]]
##
## done
##
  sumSt. <- cbind(sumSt, rmsDevOverlap=rmsd)
  attr(Dat1, 'headers') <- attr(Data[[nFiles]], 'headers')
  attr(Dat1, 'footers') <- attr(Data[[nFiles]], 'footers')
  attr(Dat1, 'summary') <- sumSt.
  Dat1
}

