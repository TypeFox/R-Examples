pValueFun <-
function(ndx, coeffMissingAllowed = 0.75) {
    notMissing <- apply(!is.na(geData[, ndx]), 1, sum)
    notMissing <- notMissing > floor((length(ndx)-1)^coeffMissingAllowed)
    aClassify <- classify(geData[notMissing, ndx])$clusters
    ans <- 1
    if(min(table(aClassify)) > floor(0.1 * nrow(geData)))
      ans <- 1 - pchisq(survdiff(stData[notMissing] ~ aClassify)$chisq, df = 1)
    return(ans)
  }
