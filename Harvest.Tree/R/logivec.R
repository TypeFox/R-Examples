logivec <- function(har.rule,label,nodenumb,newsim)
{

  r.name <- rownames(har.rule)
  logvec <- rep(TRUE, nrow(newsim))

  for (i in 1:length(r.name)){
    if(!is.infinite(har.rule[i,1]) | !is.infinite(har.rule[i,2])){
      logvec <- newsim[, r.name[i]]>= har.rule[i,1] & newsim[, r.name[i]] < har.rule[i,2] & logvec
    }
    else if(!is.na(har.rule[i,3])) {
      logvec <- sapply(newsim[,r.name[i]],function(x) grepl(x,har.rule[i,3])) & logvec
    }
    else logvec <- logvec
  }
  
  
  logvec[which(!newsim$rownn %in% nodenumb)] <- FALSE
  logvec[which(is.na(logvec))] <- newsim$rownn[which(is.na(logvec))]==label
  return(logvec)
}
