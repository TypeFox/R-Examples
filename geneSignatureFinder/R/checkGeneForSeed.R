checkGeneForSeed <-
function(geneLevels, evaluateBICs = TRUE, cutoff = 1.95)
  {
    aClassify <- classify(geneLevels)
    
    if((sum(aClassify$missing) == length(aClassify$missing)) |
         (length(unique(aClassify$clusters[!aClassify$missing])) == 1)) {
      ans <- c(NA, NA, NA, NA, rep(0, length(geneLevels)))
      names(ans) <- c("tValue", "pValue", "bic1", "bic2", names(geneLevels))
      return(ans)
    }
    
    aClassify$clusters <-  as.numeric(goodAndPoorClassification(aClassify$clusters))
    
    if(evaluateBICs)
      bics <- BICs(geneLevels[!aClassify$missing],
                   aClassify$clusters[!aClassify$missing],
                   cutoff = cutoff)$bics else bics <- c(NA, NA)
    tValue <- survdiff(stData[!aClassify$missing] ~ aClassify$clusters[!aClassify$missing])$chisq
    pValue <- 1 - pchisq(tValue, df = 1)
    ans <- c(tValue, pValue, bics, aClassify$clusters)
    names(ans) <- c("tValue", "pValue", "bic1", "bic2", names(geneLevels))
    
    return(ans)
  }
