## load packages RFLPtools (version >= 1.6), MKmisc
library(RFLPtools)
library(MKmisc)

## general parameters
M <- 5 # finally use M <- 1000
N <- 10 # increase N?
nrB <- 3:10

## helper functions
checkResultsGerm <- function(newName, resData){
  rownames(resData[[newName]])[1] == newName
}
dist2ref <- function(nB, newData, refData, distfun = dist, nrMissing = 0, LOD = 0){ 
  newData$Sample <- paste("New", newData$Sample, sep = "_")
  complData <- rbind(newData, refData)
  if(nrMissing == 0){
    res <- RFLPdist(complData, nrBands = nB, distfun = distfun, LOD = LOD)
  }else{
    res <- RFLPdist2(complData, nrBands = nB, distfun = distfun,
                     nrMissing = nrMissing)
  }
  res
}
checkResultsRFLP <- function(resData){
  resData <- as.matrix(resData)
  resData <- resData[,grepl("New", colnames(resData))]
  resData <- resData[!grepl("New", rownames(resData)),]
  newNam <- sapply(strsplit(colnames(resData), split = "_"), "[", 2)
  refNam <- rownames(resData)[apply(resData, 2, which.min)]
  sum(refNam == newNam)
}
checkResultsFragMatch <- function(resData, nrMissing = 0){
  Matches <- sapply(apply(resData, 2, strsplit, split = "_"), 
                    function(x, nrMissing){
                      sapply(x, 
                             function(x, nrMissing){ 
                               as.integer(x[1]) == (as.integer(x[2])-nrMissing)
                             }, 
                             nrMissing = nrMissing)
                    }, nrMissing = nrMissing)
  rowNams <- rownames(resData)
  colNams <- colnames(resData)
  nrMatch <- 0
  for(i in 1:ncol(Matches)){
    if(colNams[i] %in% rowNams[which(Matches[,i])]) nrMatch <- nrMatch + 1
  }
  nrMatch
}


###############################################################################
## 1. Only Measurement error
###############################################################################
acc.germ.joint1 <- numeric(M)
acc.germ.forward1 <- numeric(M)
acc.germ.backward1 <- numeric(M)
acc.germ.sum1 <- numeric(M)
acc.rflptools.eucl1 <- numeric(M)
acc.rflptools.cor1 <- numeric(M)
acc.rflptools.diff1 <- numeric(M)
acc.FragMatch1.5 <- numeric(M)
acc.FragMatch1.10 <- numeric(M)
acc.FragMatch1.25 <- numeric(M)

for(i in 1:M){
  print(i)
  # generate reference data
  refData1 <- simulateRFLPdata(N = N, bandCenters = seq(200, 800, by = 100),
                               nrBands = nrB, refData = TRUE)
  # add measurement error
  newData1 <- refData1
  newData1$MW <- newData1$MW + rnorm(nrow(newData1), mean = 0, sd = 5)
  
  ## check range of measurement error
  #summary(newData1$MW - refData1$MW)
  
  # apply germ
  res1.germ.joint <- germ(newData = newData1, refData = refData1)
  res1.germ.forward <- lapply(res1.germ.joint, function(x) x[order(x[,"Forward Max"]),])
  res1.germ.backward <- lapply(res1.germ.joint, function(x) x[order(x[,"Backward Max"]),])
  res1.germ.sum <- lapply(res1.germ.joint, function(x) x[order(x[,"Sum of Bands"]),])
  
  ## accuracy for GERM (in %)
  SampleNames <- unique(refData1$Sample)
  acc.germ.joint1[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res1.germ.joint))/length(SampleNames)
  acc.germ.forward1[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res1.germ.forward))/length(SampleNames)
  acc.germ.backward1[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res1.germ.backward))/length(SampleNames)
  acc.germ.sum1[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res1.germ.sum))/length(SampleNames)
  
  ## results for RFLPtools
  res1.rflptools.eucl <- lapply(nrB, dist2ref, newData = newData1, refData = refData1)
  res1.rflptools.cor <- lapply(nrB, dist2ref, newData = newData1, refData = refData1, dist = corDist)
  res1.rflptools.diff <- lapply(nrB, dist2ref, newData = newData1, refData = refData1, dist = diffDist)
  
  ## accuracy for RFLPtools
  acc.rflptools.eucl1[i] <- 100*sum(sapply(res1.rflptools.eucl, checkResultsRFLP))/length(SampleNames)
  acc.rflptools.cor1[i] <- 100*sum(sapply(res1.rflptools.cor, checkResultsRFLP))/length(SampleNames)
  acc.rflptools.diff1[i] <- 100*sum(sapply(res1.rflptools.diff, checkResultsRFLP))/length(SampleNames)
  
  ## apply FragMatch
  res1.FragMatch.5 <- FragMatch(newData = newData1, refData = refData1, errorBound = 5)
  res1.FragMatch.10 <- FragMatch(newData = newData1, refData = refData1, errorBound = 10)
  res1.FragMatch.25 <- FragMatch(newData = newData1, refData = refData1, errorBound = 25)
  
  ## accuracy for fraqMatch
  acc.FragMatch1.5[i] <- 100*checkResultsFragMatch(res1.FragMatch.5)/length(SampleNames)
  acc.FragMatch1.10[i] <- 100*checkResultsFragMatch(res1.FragMatch.10)/length(SampleNames)
  acc.FragMatch1.25[i] <- 100*checkResultsFragMatch(res1.FragMatch.25)/length(SampleNames)
}
mean(acc.germ.joint1)
mean(acc.germ.forward1)
mean(acc.germ.backward1)
mean(acc.germ.sum1)
mean(acc.rflptools.eucl1)
mean(acc.rflptools.cor1)
mean(acc.rflptools.diff1)
mean(acc.FragMatch1.5)
mean(acc.FragMatch1.10)
mean(acc.FragMatch1.25)


###############################################################################
## 2. Measurement error + 1 band below LOD
###############################################################################
acc.germ.joint2 <- numeric(M)
acc.germ.forward2 <- numeric(M)
acc.germ.backward2 <- numeric(M)
acc.germ.sum2 <- numeric(M)
acc.rflptools.eucl2 <- numeric(M)
acc.rflptools.cor2 <- numeric(M)
acc.rflptools.diff2 <- numeric(M)
acc.rflptools.eucl21 <- numeric(M)
acc.rflptools.cor21 <- numeric(M)
acc.rflptools.diff21 <- numeric(M)
acc.FragMatch2.5 <- numeric(M)
acc.FragMatch2.10 <- numeric(M)
acc.FragMatch2.25 <- numeric(M)

for(i in 1:M){
  print(i)
  # generate reference data
  refData2 <- simulateRFLPdata(N = N, bandCenters = seq(200, 800, by = 100),
                               nrBands = nrB, refData = TRUE)
  # add band below LOD
  SampleNames <- unique(refData2$Sample)
  for(j in 1:length(SampleNames)){
    temp <- refData2[refData2$Sample == SampleNames[j],]
    addLOD <- data.frame(Sample = SampleNames[j],
                         Band = 0,
                         MW = runif(1, min = 0, max = 100),
                         Enzyme = temp$Enzyme[1],
                         Taxonname = temp$Taxonname[1],
                         Accession = temp$Accession[1])
    temp <- rbind(addLOD, temp)
    temp$Band <- temp$Band + 1
    if(j == 1) 
      newData2 <- temp
    else
      newData2 <- rbind(newData2, temp)
  }
  rownames(newData2) <- 1:nrow(newData2)
  
  # add measurement error
  newData2$MW <- newData2$MW + rnorm(nrow(newData2), mean = 0, sd = 5)
  
  # apply germ
  res2.germ.joint <- germ(newData = newData2, refData = refData2)
  res2.germ.forward <- lapply(res2.germ.joint, function(x) x[order(x[,"Forward Max"]),])
  res2.germ.backward <- lapply(res2.germ.joint, function(x) x[order(x[,"Backward Max"]),])
  res2.germ.sum <- lapply(res2.germ.joint, function(x) x[order(x[,"Sum of Bands"]),])
  
  ## accuracy for GERM (in %)
  acc.germ.joint2[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res2.germ.joint))/length(SampleNames)
  acc.germ.forward2[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res2.germ.forward))/length(SampleNames)
  acc.germ.backward2[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res2.germ.backward))/length(SampleNames)
  acc.germ.sum2[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res2.germ.sum))/length(SampleNames)
  
  ## results for RFLPtools
  res2.rflptools.eucl <- lapply(nrB, dist2ref, newData = newData2, refData = refData2, LOD = 125)
  res2.rflptools.cor <- lapply(nrB, dist2ref, newData = newData2, refData = refData2, dist = corDist, LOD = 125)
  res2.rflptools.diff <- lapply(nrB, dist2ref, newData = newData2, refData = refData2, dist = diffDist, LOD = 125)
  res2.rflptools.eucl1 <- lapply(nrB, dist2ref, newData = newData2, refData = refData2, nrMissing = 1)
  res2.rflptools.cor1 <- lapply(nrB, dist2ref, newData = newData2, refData = refData2, dist = corDist, nrMissing = 1)
  res2.rflptools.diff1 <- lapply(nrB, dist2ref, newData = newData2, refData = refData2, dist = diffDist, nrMissing = 1)
  
  ## accuracy for RFLPtools
  acc.rflptools.eucl2[i] <- 100*sum(sapply(res2.rflptools.eucl, checkResultsRFLP))/length(SampleNames)
  acc.rflptools.cor2[i] <- 100*sum(sapply(res2.rflptools.cor, checkResultsRFLP))/length(SampleNames)
  acc.rflptools.diff2[i] <- 100*sum(sapply(res2.rflptools.diff, checkResultsRFLP))/length(SampleNames)
  acc.rflptools.eucl21[i] <- 100*sum(sapply(res2.rflptools.eucl1, checkResultsRFLP))/length(SampleNames)
  acc.rflptools.cor21[i] <- 100*sum(sapply(res2.rflptools.cor1, checkResultsRFLP))/length(SampleNames)
  acc.rflptools.diff21[i] <- 100*sum(sapply(res2.rflptools.diff1, checkResultsRFLP))/length(SampleNames)
  
  ## apply FragMatch
  res2.FragMatch.5 <- FragMatch(newData = newData2, refData = refData2, errorBound = 5)
  res2.FragMatch.10 <- FragMatch(newData = newData2, refData = refData2, errorBound = 10)
  res2.FragMatch.25 <- FragMatch(newData = newData2, refData = refData2, errorBound = 25)
  
  ## accuracy for fraqMatch
  acc.FragMatch2.5[i] <- 100*checkResultsFragMatch(res2.FragMatch.5)/length(SampleNames)
  acc.FragMatch2.10[i] <- 100*checkResultsFragMatch(res2.FragMatch.10)/length(SampleNames)
  acc.FragMatch2.25[i] <- 100*checkResultsFragMatch(res2.FragMatch.25)/length(SampleNames)
}
mean(acc.germ.joint2)
mean(acc.germ.forward2)
mean(acc.germ.backward2)
mean(acc.germ.sum2)
mean(acc.rflptools.eucl2)
mean(acc.rflptools.cor2)
mean(acc.rflptools.diff2)
mean(acc.rflptools.eucl21)
mean(acc.rflptools.cor21)
mean(acc.rflptools.diff21)
mean(acc.FragMatch2.5)
mean(acc.FragMatch2.10)
mean(acc.FragMatch2.25)


###############################################################################
## 3. Measurement error + 1 faint band
###############################################################################
acc.germ.joint3 <- numeric(M)
acc.germ.forward3 <- numeric(M)
acc.germ.backward3 <- numeric(M)
acc.germ.sum3 <- numeric(M)
acc.rflptools.eucl3 <- numeric(M)
acc.rflptools.cor3 <- numeric(M)
acc.rflptools.diff3 <- numeric(M)
acc.FragMatch3.5 <- numeric(M)
acc.FragMatch3.10 <- numeric(M)
acc.FragMatch3.25 <- numeric(M)

for(i in 1:M){
  print(i)
  # generate reference data
  refData3 <- simulateRFLPdata(N = N, bandCenters = seq(200, 800, by = 100),
                               nrBands = nrB, refData = TRUE)
  # add one addtional band "somewhere"
  SampleNames <- unique(refData3$Sample)
  for(j in 1:length(SampleNames)){
    temp <- refData3[refData3$Sample == SampleNames[j],]
    addBand <- data.frame(Sample = SampleNames[j],
                         Band = max(temp$Band)+1,
                         MW = -runif(1, min = 150, max = 850),
                         Enzyme = temp$Enzyme[1],
                         Taxonname = temp$Taxonname[1],
                         Accession = temp$Accession[1])
    temp <- rbind(temp, addBand)
    temp <- temp[order(abs(temp$MW)),]
    temp$Band <- 1:nrow(temp)
    if(j == 1) 
      newData3 <- temp
    else
      newData3 <- rbind(newData3, temp)
  }
  rownames(newData3) <- 1:nrow(newData3)
  
  # add measurement error
  newData3$MW <- newData3$MW + rnorm(nrow(newData3), mean = 0, sd = 5)
  
  # apply germ
  res3.germ.joint <- germ(newData = newData3, refData = refData3)
  res3.germ.forward <- lapply(res3.germ.joint, function(x) x[order(x[,"Forward Max"]),])
  res3.germ.backward <- lapply(res3.germ.joint, function(x) x[order(x[,"Backward Max"]),])
  res3.germ.sum <- lapply(res3.germ.joint, function(x) x[order(x[,"Sum of Bands"]),])
  
  ## accuracy for GERM (in %)
  acc.germ.joint3[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res3.germ.joint))/length(SampleNames)
  acc.germ.forward3[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res3.germ.forward))/length(SampleNames)
  acc.germ.backward3[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res3.germ.backward))/length(SampleNames)
  acc.germ.sum3[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res3.germ.sum))/length(SampleNames)
  
  ## results for RFLPtools
  newData31 <- newData3
  newData31$MW <- abs(newData31$MW)
  res3.rflptools.eucl1 <- lapply(nrB, dist2ref, newData = newData31, refData = refData3, nrMissing = 1)
  res3.rflptools.cor1 <- lapply(nrB, dist2ref, newData = newData31, refData = refData3, dist = corDist, nrMissing = 1)
  res3.rflptools.diff1 <- lapply(nrB, dist2ref, newData = newData31, refData = refData3, dist = diffDist, nrMissing = 1)
  
  ## accuracy for RFLPtools
  acc.rflptools.eucl3[i] <- 100*sum(sapply(res3.rflptools.eucl1, checkResultsRFLP))/length(SampleNames)
  acc.rflptools.cor3[i] <- 100*sum(sapply(res3.rflptools.cor1, checkResultsRFLP))/length(SampleNames)
  acc.rflptools.diff3[i] <- 100*sum(sapply(res3.rflptools.diff1, checkResultsRFLP))/length(SampleNames)
  
  ## apply FragMatch
  res3.FragMatch.5 <- FragMatch(newData = newData31, refData = refData3, errorBound = 5)
  res3.FragMatch.10 <- FragMatch(newData = newData31, refData = refData3, errorBound = 10)
  res3.FragMatch.25 <- FragMatch(newData = newData31, refData = refData3, errorBound = 25)
  
  ## accuracy for fraqMatch
  acc.FragMatch3.5[i] <- 100*checkResultsFragMatch(res3.FragMatch.5)/length(SampleNames)
  acc.FragMatch3.10[i] <- 100*checkResultsFragMatch(res3.FragMatch.10)/length(SampleNames)
  acc.FragMatch3.25[i] <- 100*checkResultsFragMatch(res3.FragMatch.25)/length(SampleNames)
}
mean(acc.germ.joint3)
mean(acc.germ.forward3)
mean(acc.germ.backward3)
mean(acc.germ.sum3)
mean(acc.rflptools.eucl3)
mean(acc.rflptools.cor3)
mean(acc.rflptools.diff3)
mean(acc.FragMatch3.5)
mean(acc.FragMatch3.10)
mean(acc.FragMatch3.25)


###############################################################################
## 4. Measurement error + 1 missing band
###############################################################################
acc.germ.joint4 <- numeric(M)
acc.germ.forward4 <- numeric(M)
acc.germ.backward4 <- numeric(M)
acc.germ.sum4 <- numeric(M)
acc.rflptools.eucl4 <- numeric(M)
acc.rflptools.cor4 <- numeric(M)
acc.rflptools.diff4 <- numeric(M)
acc.FragMatch4.5 <- numeric(M)
acc.FragMatch4.10 <- numeric(M)
acc.FragMatch4.25 <- numeric(M)

for(i in 1:M){
  print(i)
  # generate reference data
  refData4 <- simulateRFLPdata(N = N, bandCenters = seq(200, 800, by = 100),
                               nrBands = nrB, refData = TRUE)
  # remove one band randomly
  SampleNames <- unique(refData4$Sample)
  for(j in 1:length(SampleNames)){
    temp <- refData4[refData4$Sample == SampleNames[j],]
    temp <- temp[-sample(1:nrow(temp), 1),]
    temp$Band <- 1:nrow(temp)
    if(j == 1) 
      newData4 <- temp
    else
      newData4 <- rbind(newData4, temp)
  }
  rownames(newData4) <- 1:nrow(newData4)
  
  # add measurement error
  newData4$MW <- newData4$MW + rnorm(nrow(newData4), mean = 0, sd = 5)
  
  # apply germ
  res4.germ.joint <- germ(newData = newData4, refData = refData4)
  res4.germ.forward <- lapply(res4.germ.joint, function(x) x[order(x[,"Forward Max"]),])
  res4.germ.backward <- lapply(res4.germ.joint, function(x) x[order(x[,"Backward Max"]),])
  res4.germ.sum <- lapply(res4.germ.joint, function(x) x[order(x[,"Sum of Bands"]),])
  
  ## accuracy for GERM (in %)
  acc.germ.joint4[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res4.germ.joint))/length(SampleNames)
  acc.germ.forward4[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res4.germ.forward))/length(SampleNames)
  acc.germ.backward4[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res4.germ.backward))/length(SampleNames)
  acc.germ.sum4[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res4.germ.sum))/length(SampleNames)
  
  ## results for RFLPtools
  res4.rflptools.eucl1 <- lapply(nrB, dist2ref, newData = newData4, refData = refData4, nrMissing = 1)
  res4.rflptools.cor1 <- lapply(nrB, dist2ref, newData = newData4, refData = refData4, dist = corDist, nrMissing = 1)
  res4.rflptools.diff1 <- lapply(nrB, dist2ref, newData = newData4, refData = refData4, dist = diffDist, nrMissing = 1)
  
  ## accuracy for RFLPtools
  acc.rflptools.eucl4[i] <- 100*sum(sapply(res4.rflptools.eucl1[-length(nrB)], checkResultsRFLP))/length(SampleNames)
  acc.rflptools.cor4[i] <- 100*sum(sapply(res4.rflptools.cor1[-length(nrB)], checkResultsRFLP))/length(SampleNames)
  acc.rflptools.diff4[i] <- 100*sum(sapply(res4.rflptools.diff1[-length(nrB)], checkResultsRFLP))/length(SampleNames)
  
  ## apply FragMatch
  res4.FragMatch.5 <- FragMatch(newData = newData4, refData = refData4, errorBound = 5)
  res4.FragMatch.10 <- FragMatch(newData = newData4, refData = refData4, errorBound = 10)
  res4.FragMatch.25 <- FragMatch(newData = newData4, refData = refData4, errorBound = 25)
  
  ## accuracy for fraqMatch
  acc.FragMatch4.5[i] <- 100*checkResultsFragMatch(res4.FragMatch.5, nrMissing = 1)/length(SampleNames)
  acc.FragMatch4.10[i] <- 100*checkResultsFragMatch(res4.FragMatch.10, nrMissing = 1)/length(SampleNames)
  acc.FragMatch4.25[i] <- 100*checkResultsFragMatch(res4.FragMatch.25, nrMissing = 1)/length(SampleNames)
}
mean(acc.germ.joint4)
mean(acc.germ.forward4)
mean(acc.germ.backward4)
mean(acc.germ.sum4)
mean(acc.rflptools.eucl4)
mean(acc.rflptools.cor4)
mean(acc.rflptools.diff4)
mean(acc.FragMatch4.5)
mean(acc.FragMatch4.10)
mean(acc.FragMatch4.25)


###############################################################################
## 5. Measurement error + 1 replicated band
###############################################################################
acc.germ.joint5 <- numeric(M)
acc.germ.forward5 <- numeric(M)
acc.germ.backward5 <- numeric(M)
acc.germ.sum5 <- numeric(M)
acc.rflptools.eucl5 <- numeric(M)
acc.rflptools.cor5 <- numeric(M)
acc.rflptools.diff5 <- numeric(M)
acc.FragMatch5.5 <- numeric(M)
acc.FragMatch5.10 <- numeric(M)
acc.FragMatch5.25 <- numeric(M)

for(i in 1:M){
  print(i)
  # generate reference data
  refData5 <- simulateRFLPdata(N = N, bandCenters = seq(200, 800, by = 100),
                               nrBands = nrB, refData = TRUE)
  # select randomly one band, replicate this band, add measurement error
  SampleNames <- unique(refData5$Sample)
  for(j in 1:length(SampleNames)){
    temp <- refData5[refData5$Sample == SampleNames[j],]
    addBand <- temp[sample(1:nrow(temp), 1),, drop = FALSE]
    temp <- rbind(temp, addBand)
    temp <- temp[order(abs(temp$MW)),]
    temp$Band <- 1:nrow(temp)
    if(j == 1) 
      newData5 <- temp
    else
      newData5 <- rbind(newData5, temp)
  }
  rownames(newData5) <- 1:nrow(newData5)
  
  # add measurement error
  newData5$MW <- newData5$MW + rnorm(nrow(newData5), mean = 0, sd = 5)
  
  # apply germ
  res5.germ.joint <- germ(newData = newData5, refData = refData5)
  res5.germ.forward <- lapply(res5.germ.joint, function(x) x[order(x[,"Forward Max"]),])
  res5.germ.backward <- lapply(res5.germ.joint, function(x) x[order(x[,"Backward Max"]),])
  res5.germ.sum <- lapply(res5.germ.joint, function(x) x[order(x[,"Sum of Bands"]),])
  
  ## accuracy for GERM (in %)
  acc.germ.joint5[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res5.germ.joint))/length(SampleNames)
  acc.germ.forward5[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res5.germ.forward))/length(SampleNames)
  acc.germ.backward5[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res5.germ.backward))/length(SampleNames)
  acc.germ.sum5[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res5.germ.sum))/length(SampleNames)
  
  ## results for RFLPtools
  res5.rflptools.eucl1 <- lapply(nrB, dist2ref, newData = newData5, refData = refData5, nrMissing = 1)
  res5.rflptools.cor1 <- lapply(nrB, dist2ref, newData = newData5, refData = refData5, dist = corDist, nrMissing = 1)
  res5.rflptools.diff1 <- lapply(nrB, dist2ref, newData = newData5, refData = refData5, dist = diffDist, nrMissing = 1)
  
  ## accuracy for RFLPtools
  acc.rflptools.eucl5[i] <- 100*sum(sapply(res5.rflptools.eucl1, checkResultsRFLP))/length(SampleNames)
  acc.rflptools.cor5[i] <- 100*sum(sapply(res5.rflptools.cor1, checkResultsRFLP))/length(SampleNames)
  acc.rflptools.diff5[i] <- 100*sum(sapply(res5.rflptools.diff1, checkResultsRFLP))/length(SampleNames)
  
  ## apply FragMatch
  res5.FragMatch.5 <- FragMatch(newData = newData5, refData = refData5, errorBound = 5)
  res5.FragMatch.10 <- FragMatch(newData = newData5, refData = refData5, errorBound = 10)
  res5.FragMatch.25 <- FragMatch(newData = newData5, refData = refData5, errorBound = 25)
  
  ## accuracy for fraqMatch
  acc.FragMatch5.5[i] <- 100*checkResultsFragMatch(res5.FragMatch.5)/length(SampleNames)
  acc.FragMatch5.10[i] <- 100*checkResultsFragMatch(res5.FragMatch.10)/length(SampleNames)
  acc.FragMatch5.25[i] <- 100*checkResultsFragMatch(res5.FragMatch.25)/length(SampleNames)
}
mean(acc.germ.joint5)
mean(acc.germ.forward5)
mean(acc.germ.backward5)
mean(acc.germ.sum5)
mean(acc.rflptools.eucl5)
mean(acc.rflptools.cor5)
mean(acc.rflptools.diff5)
mean(acc.FragMatch5.5)
mean(acc.FragMatch5.10)
mean(acc.FragMatch5.25)


###############################################################################
## 6. Measurement error + positive shift of bands
###############################################################################
# small shifts are no problem for the algorithms -> use larger shift
shift <- 50
acc.germ.joint6 <- numeric(M)
acc.germ.forward6 <- numeric(M)
acc.germ.backward6 <- numeric(M)
acc.germ.sum6 <- numeric(M)
acc.rflptools.eucl6 <- numeric(M)
acc.rflptools.cor6 <- numeric(M)
acc.rflptools.diff6 <- numeric(M)
acc.FragMatch6.5 <- numeric(M)
acc.FragMatch6.10 <- numeric(M)
acc.FragMatch6.25 <- numeric(M)

for(i in 1:M){
  print(i)
  # generate reference data
  refData6 <- simulateRFLPdata(N = N, bandCenters = seq(200, 800, by = 100),
                               nrBands = nrB, refData = TRUE)
  # fixed shift of bands
  newData6 <- refData6
  newData6$MW <- newData6$MW + shift
  
  # add measurement error
  newData6$MW <- newData6$MW + rnorm(nrow(newData6), mean = 0, sd = 5)
  
  ## check range of measurement error
  #summary(newData6$MW - refData6$MW)
  
  # apply germ
  res6.germ.joint <- germ(newData = newData6, refData = refData6)
  res6.germ.forward <- lapply(res6.germ.joint, function(x) x[order(x[,"Forward Max"]),])
  res6.germ.backward <- lapply(res6.germ.joint, function(x) x[order(x[,"Backward Max"]),])
  res6.germ.sum <- lapply(res6.germ.joint, function(x) x[order(x[,"Sum of Bands"]),])
  
  ## accuracy for GERM (in %)
  SampleNames <- unique(refData6$Sample)
  acc.germ.joint6[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res6.germ.joint))/length(SampleNames)
  acc.germ.forward6[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res6.germ.forward))/length(SampleNames)
  acc.germ.backward6[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res6.germ.backward))/length(SampleNames)
  acc.germ.sum6[i] <- 100*sum(sapply(SampleNames, checkResultsGerm, resData = res6.germ.sum))/length(SampleNames)
  
  ## results for RFLPtools
  res6.rflptools.eucl <- lapply(nrB, dist2ref, newData = newData6, refData = refData6)
  res6.rflptools.cor <- lapply(nrB, dist2ref, newData = newData6, refData = refData6, dist = corDist)
  res6.rflptools.diff <- lapply(nrB, dist2ref, newData = newData6, refData = refData6, dist = diffDist)
  
  ## accuracy for RFLPtools
  acc.rflptools.eucl6[i] <- 100*sum(sapply(res6.rflptools.eucl, checkResultsRFLP))/length(SampleNames)
  acc.rflptools.cor6[i] <- 100*sum(sapply(res6.rflptools.cor, checkResultsRFLP))/length(SampleNames)
  acc.rflptools.diff6[i] <- 100*sum(sapply(res6.rflptools.diff, checkResultsRFLP))/length(SampleNames)
  
  ## apply FragMatch
  res6.FragMatch.5 <- FragMatch(newData = newData6, refData = refData6, errorBound = 5)
  res6.FragMatch.10 <- FragMatch(newData = newData6, refData = refData6, errorBound = 10)
  res6.FragMatch.25 <- FragMatch(newData = newData6, refData = refData6, errorBound = 25)
  
  ## accuracy for fraqMatch
  acc.FragMatch6.5[i] <- 100*checkResultsFragMatch(res6.FragMatch.5)/length(SampleNames)
  acc.FragMatch6.10[i] <- 100*checkResultsFragMatch(res6.FragMatch.10)/length(SampleNames)
  acc.FragMatch6.25[i] <- 100*checkResultsFragMatch(res6.FragMatch.25)/length(SampleNames)
}
mean(acc.germ.joint6)
mean(acc.germ.forward6)
mean(acc.germ.backward6)
mean(acc.germ.sum6)
mean(acc.rflptools.eucl6)
mean(acc.rflptools.cor6)
mean(acc.rflptools.diff6)
mean(acc.FragMatch6.5)
mean(acc.FragMatch6.10)
mean(acc.FragMatch6.25)

