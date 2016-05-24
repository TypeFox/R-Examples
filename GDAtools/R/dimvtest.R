dimvtest <- function(resmca,l,n=names(l),dim=1:resmca$call$ncp) {
  res1 <- list()
  res2 <- list()
  for(i in 1:length(l)) {
    res1[[i]] <- varsup(resmca,l[[i]])$v.test
    rownames(res1[[i]]) <- paste(n[[i]],rownames(res1[[i]]),sep='.')
    res2[[i]] <- varsup(resmca,l[[i]])$weight
    }
  res1 <- do.call('rbind.data.frame',res1)
  res2 <- unlist(res2)
  res <- list()
  for(i in 1:length(dim)) {
    z <- data.frame(res1[,i],res2)
    colnames(z) <- c('vtest','weight')
    rownames(z) <- rownames(res1)
    z <- z[order(-z$vtest),]
    z <- z[abs(z$vtest)>=2.58,]
    res[[i]] <- z
    }
  names(res) <- paste('dim',dim,sep='.')
  return(res)
  }