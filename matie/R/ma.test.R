ma.test <-
function (d, maStruct, permPartition=maStruct$partition, lookup=TRUE, reps=1000)
{
  
  pal <-
    function (maStruct)
    {
      
      if(maStruct$mCols != 2) stop("lookup only valid for 2D data sets")
      if(maStruct$nRows > 400) warning(paste("nRows =",toString(maStruct$nRows),
                                             ", lookup only robust for nRows < 400."))
      
      LRTsampleSizes <- c(12,25,50,100,200,400)
      
      LRTweights <- rbind(c(0.581204, 0.33861, 0.011954, 0.0682328),
                          c(0.552741, 0.0735269, 0.373732, 0.0),
                          c(0.555714, 0.0127543, 0.373861, 0.0576714),
                          c(0.595711, 0.346582, 0.0206349, 0.0370724),
                          c(0.663338, 0.308716, 0., 0.027946),
                          c(0.71685, 0.262288, 7.951*10^-9, 0.0208622))
      
      LRTdfs <- rbind(c(1.46232, 1.65733*10^-6, 4.08098),
                      c(4.9455, 1.37708, 2.81088),
                      c(7.69287, 1.5031, 4.67881),
                      c(2.0419, 0.0000589314, 7.17571),
                      c(2.47259, 0.267051, 8.12206),
                      c(2.94421, 0.00177885, 9.1659))
      
      knownParamP <- function (LRS, weights, dfs)
      {
        pVals = sum((1-pchisq(LRS, dfs))*weights[2:length(weights)])
        return(pVals)
      }
      
      pValue <- function (LRS, n)
      {
        pValVector <- LRTsampleSizes*0
        for (i in 1:length(pValVector))
        {
          pValVector[i] <- knownParamP(LRS,LRTweights[i,],LRTdfs[i,])
        }
        logSampleSizes <- log(LRTsampleSizes)
        p<-approx(logSampleSizes,pValVector,log(n),rule = 2:2)$y
        return(p)
      }
      
      return(pValue(maStruct$LRstat,maStruct$nRows))
      
    }
  
  if(lookup) return(pal(maStruct))
  
  # no lookup, so we use reps to perform monte carlo simulation
  
  testPartition <- maStruct$partition
  
  permData <- function (dat,partition)
  {
    rows <- nrow(dat)
    return(do.call(cbind,lapply(partition,(function(x)cbind(dat[,x])[sample(1:rows),]))))
  }
  
  permECDF <- function (dat,permpartition,testpartition,reps)
  {
    return(ecdf(replicate(reps,ma(permData(dat,permpartition),partition=testpartition)$LRstat)))
  }
  
  
  LRecdf <- permECDF(d,permPartition,testPartition,reps=reps)
  datLRS <- maStruct$LRstat
  return(1-LRecdf(datLRS))
  
}
