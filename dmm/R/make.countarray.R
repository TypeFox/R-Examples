make.countarray <-
function(mdf,traits){
# make.countarray() - do countnotna() on all traitpairs and store counts in matrix returned
  if(is.null(mdf$rel)) {
    df <- mdf
  } else {
    df <- mdf$df
  }

  l <- length(traits)
  counts <- array(0,c(l,l))
  dimnames(counts) <- list(traits,traits)
  for(i in traits){
    counti <- countnotna(df[,i])
    for(j in traits){
      countj <- countnotna(df[,j])
      counts[i,j] <- length(df[,i][!is.na(df[,i]) & !is.na(df[,j])])
    }
  }
  return(counts)
}
