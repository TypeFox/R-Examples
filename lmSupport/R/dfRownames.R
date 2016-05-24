dfRownames <-
function(Data, SubID = 'SubID', FixedWidth = TRUE, Remove=TRUE, MaxNumDigits = NULL)
#Sets the row names of the data frame to the variable name listed as SubID.
#SubID should be text name of variable.  
#Also keeps number of characters constant by default and removes SubID by default
{
  TheRowNames = Data[,SubID]
  
  #determine number of digits for SubID.  If null, set to max observed length
  if(is.null(MaxNumDigits)) MaxNumDigits = floor(log10(max(TheRowNames)))+1 
  
  if (FixedWidth)  
    {
       FixedNames= NA
       for (i in 1:nrow(Data))  FixedNames[i] = paste(paste(rep(0,MaxNumDigits-(floor(log10(TheRowNames[i]))+1)),collapse=''),TheRowNames[i], sep = '')
       TheRowNames = FixedNames
    }
  rownames(Data) = TheRowNames
  if (Remove) Data[,SubID] = NULL 
  return(Data)
}

