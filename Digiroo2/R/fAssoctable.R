fAssoctable <-
function(dnndata)
{
  require(spdep)
  
  animID <- attr(dnndata, "region.id")  
  idata4 <- fextractassoc(dnndata)
  if(all(idata4[[1]][1]==0)) idata4 <- idata4[-1] #Removes 0 value if present in list object cell 1
  #idata4 <- lapply(lapply(new.assoc4, function(x){replace(x, x == 0, NA)}), function(x) x[!is.na(x)])
  list.length <- order(unlist(lapply(idata4,length)),decreasing = TRUE)
  idata5 <- lapply(1:length(list.length),function(i) idata4[[list.length[i]]])
  
  df1 <- data.frame(Group=0,IDs=0,stringsAsFactors=FALSE)[NULL,]
  iCounter <- 1
  usedIDs <- 0
  
  while(length(unique(usedIDs)) < length(animID) & iCounter <= length(idata5))
  {
    (IDs <- idata5[[iCounter]])
    if(iCounter > 1)
    {
      if(length(intersect(IDs, usedIDs)) == 0)
      {
        Group <-  Group + 1
        df2 <- cbind(Group,IDs)
        df1 <- rbind(df1,df2)
        usedIDs <- c(usedIDs,IDs)
      }
      if(length(setdiff(IDs,usedIDs)) > 0)
      {
        df2 <- cbind(Group,setdiff(IDs,usedIDs))
        df1 <- rbind(df1,df2)    
        usedIDs <- c(usedIDs,IDs)
      }
    }else{
      Group <- 1
      df1 <- cbind(Group,IDs)
      usedIDs <- IDs
    }
    iCounter <- iCounter + 1
  }
  df1[,2] <- as.character(animID[df1[,2]])  #Change the IDs from factor level numbers back into ID names
  return(df1)
}
