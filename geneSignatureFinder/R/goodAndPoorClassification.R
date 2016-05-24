goodAndPoorClassification <-
function (clustering) 
{
  if (areDataNotLoaded()) 
    return(NULL)
  
  if(!is.factor(clustering))
	  clustering <- as.factor(clustering)
  
  missing <- is.na(clustering)
  x <- as.numeric(clustering[!missing])
  levs <- sort(unique(x))
  tmp <- survdiff(stData[!missing] ~ x)
  
  ans <- rep(NA, length(clustering))
  ans[!missing] <- "poor"
  if(tmp$obs[1] < tmp$exp[1]) 
    ans[!missing][x == levs[1]] <- "good" else 
      ans[!missing][x == levs[2]] <- "good"
  return(as.factor(ans))
}
