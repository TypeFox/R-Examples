gSquare = function(target, dataset, xIndex, csIndex, dataInfo=NULL, univariateModels=NULL, hash = FALSE, 
stat_hash=NULL, pvalue_hash=NULL,robust=FALSE) {
  #Conditional Independence test based on the G test of independence (log likelihood ratio  test)
  
  csIndex[which(is.na(csIndex))] = 0;
  
  if(hash == TRUE)
  {
    csIndex2 = csIndex[which(csIndex!=0)]
    csindex2 = sort(csIndex2)
    xcs = c(xIndex,csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if(is.null(stat_hash[[key]]) == FALSE)
    {
      stat = stat_hash[[key]];
      pvalue = pvalue_hash[[key]];
      flag = 1;
      
      results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  
  xIndex = as.integer(xIndex);
  csIndex = as.integer(csIndex);
  if(length(csIndex) == 1)
  {
    if(xIndex == csIndex)
    {
      if(hash == TRUE)#update hash objects
      {
        stat_hash[[key]] <- 0;#.set(stat_hash , key , 0)
        pvalue_hash[[key]] <- log(1);#.set(pvalue_hash , key , 1)
      }
      pvalue = log(1);
      stat = 0;
      flag = 1;
      results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
    
  }else
  {
    csIndex = csIndex[csIndex!=0]
    if(length(csIndex) == 0)
    {
      csIndex = NULL
    }
  }
  
  if(is.null(csIndex))
  {
    #if the univariate models have been already compute
    if(!is.null(univariateModels))
    {
      pvalue = univariateModels$pvalue[[xIndex]];
      stat = univariateModels$stat[[xIndex]];
      flag = univariateModels$flag[[xIndex]];
      results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  p <- ncol(dataset) + 1
  #levels for each variable
   mod <- cat.ci(p, xIndex, csIndex, cbind(dataset, target) )
     stat <- mod[1]
     pvalue <- mod[2]
  flag = 1;
  #define a temp statistic due to the pvalue because gSquareBin does not return the stat

  
  if(hash == TRUE)
  {
    stat_hash[[key]] <- stat;#.set(stat_hash , key , stat)
    pvalue_hash[[key]] <- pvalue;#.set(pvalue_hash , key , pvalue)
  }
  
  results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
}








