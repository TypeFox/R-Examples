tValueFun <-
function(ndx, coeffMissingAllowed = 0.75) {
    notMissing <- apply(!is.na(geData[, ndx]), 1, sum)
    notMissing <- notMissing > floor((length(ndx)-1)^coeffMissingAllowed)
    clusters <- classify(geData[notMissing, ndx])$clusters
    ans <- 0 #variato
    if(min(table(clusters)) > floor(0.1 * nrow(geData)))
      ans <- survdiff(stData[notMissing] ~ clusters)$chisq #variato
    return(ans)
  }
