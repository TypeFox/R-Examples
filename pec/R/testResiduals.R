testResiduals <- function(object,
                          times,
                          testTimes,
                          rangeInt,
                          confInt,
                          confLevel,
                          keepTestedResiduals){
  NF <- length(object)
  comparisonList <- allComparisons(names(object))
  testExact <- NROW(object[[1]])<100
  # {{{ compute integrated residuals 
  
  testIBS <- !is.null(rangeInt)
  if (testIBS==TRUE){
    if (length(rangeInt)==2 && is.numeric(rangeInt))
      range <- rangeInt
    else
      range <- NULL
    integratedResiduals <- lapply(object,function(x){
      apply(x,1,function(r){
        Dint(x=times,y=r,range=range,restrictNonMissing=FALSE)
      })})
    ## naFractionIBS <- lapply(integratedResiduals,function(x)mean(is.na(x)))
  }
  # }}}
  # {{{ extract residuals at testTimes
  if (!is.null(testTimes)){
    timePos <- prodlim::sindex(times,testTimes)
    testTimeResiduals <- lapply(object,function(x){
      x[,timePos,drop=FALSE]
    })
    ## naFractionTestTimes <- lapply(testTimeResiduals,function(x)colMeans(is.na(x)))
  }
  # }}}
  loop <- lapply(comparisonList,function(cc){
    # {{{ test residuals at time points
    if (!is.null(testTimes)){
      Rdiff <- testTimeResiduals[[cc[2]]]-testTimeResiduals[[cc[1]]]
      wtest <- lapply(1:length(testTimes),function(t){
        d <- Rdiff[,t,drop=TRUE]
        if (any(is.na(d))){
          list(p.value=NA,conf.int=c(NA,NA))
        }
        else{
          suppressWarnings(wilcox.test(d,alternative="less",exact=testExact,conf.int=confInt,conf.level=confLevel))
        }
      })
      loopOut <- list(pValue=sapply(wtest,function(w)w$p.value))
      if (confInt==TRUE){
        loopOut <- c(loopOut,list(upperLimit=sapply(wtest,function(w)w$conf.int[2])))
      }
    }
    else{
      loopOut <- vector(mode = "list", length = NF)
    }
    # }}}
    # {{{ test integrated residuals
    if (testIBS){
      dIBS <- integratedResiduals[[cc[2]]]-integratedResiduals[[cc[1]]]
      if (any(is.na(dIBS))){
        loopOut <- c(loopOut,list(IBSpValue=NA))
        if (confInt==TRUE){
          loopOut <- c(loopOut,list(IBSupper=NA))
        }
      }
      else{
        wtestIBS <- suppressWarnings(wilcox.test(dIBS,alternative="less",exact=testExact,conf.int=confInt,conf.level=confLevel))
        loopOut <- c(loopOut,list(IBSpValue=wtestIBS$p.value))
        if (confInt==TRUE){
          loopOut <- c(loopOut,list(IBSupper=wtestIBS$conf.int[2]))
        }
      }
    }
    # }}}
    loopOut
  })
  # {{{ prepare output
  if (!is.null(testTimes)){
    out <- list(pValues=lapply(loop,function(x)x$pValue))
  }
  else{
    out <- NULL
  }
  if (testIBS){
    out <- c(out,list(IBSpValue=lapply(loop,function(x)x$IBSpValue)))
  }
  # }}}
  out
}
  
