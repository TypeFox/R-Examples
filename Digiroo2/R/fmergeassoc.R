fmergeassoc <-
function(idata2)
{
  counter <- 1
  inewdata <- list(0)
  for(j in 1:length(idata2))
  {
    for (i in 1:length(idata2))
    {
      if(length(intersect(idata2[[j]],idata2[[i]])) > 0)
        if(all(is.element(idata2[[j]],idata2[[i]])) == FALSE)
        {
          inewdata[[counter]] <- union(idata2[[j]],idata2[[i]])
          counter <- counter + 1
        }
    }
  }
  return(unique(lapply(c(inewdata,idata2),sort)))
}
