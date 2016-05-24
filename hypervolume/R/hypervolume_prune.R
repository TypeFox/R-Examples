hypervolume_prune <- function(hvlist, minnp=NULL, minvol=NULL, returnids=FALSE)
{
  if (!is.null(minnp))
  {
    if (minnp < 0)
    {
      stop("minnp must be at least zero.")
    }
    if (!is.null(minvol))
    {
      stop("Cannot specify both minvol and minnp.")
    }
  }
  else if (!is.null(minvol))
  {
    if (minvol < 0)
    {
      stop("minvol must be at least zero.")
    }
    
    if (!is.null(minnp))
    {
      stop("Cannot specify both minvol and minnp.")
    }
  }
  else
  {
    stop("Must specify either minvol or minnp")
  }  
  
  if(class(hvlist) != "HypervolumeList")
  {
    stop("Input hvlist must be of class HypervolumeList.")
  }
  
  # do segmentation
  dodrop <- rep(FALSE, length(hvlist@HVList))
  for (i in 1:length(hvlist@HVList))
  {
    np <- nrow(hvlist@HVList[[i]]@RandomUniformPointsThresholded)
    vol <- hvlist@HVList[[i]]@Volume
    
    if (!is.null(minnp))
    {
      if (np < minnp)
      {
        dodrop[i] <- TRUE
      }
    }
    if (!is.null(minvol))
    {
      if (vol < minvol)
      {
        dodrop[i] <- TRUE
      }
    }
  }
  
  hvlist@HVList <- hvlist@HVList[!dodrop]
  
  if (returnids)
  {
    return(list(HVList=hvlist, IDs=which(!dodrop)))
  }
  else
  {
    return(hvlist)
  }
}