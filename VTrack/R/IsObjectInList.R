IsObjectInList <-
function(TheList, TheObject)
{
    fObjectInList <- FALSE
    
    if (length(TheList) > 0)
        for (i in 1:length(TheList))
            if (TheList[i] == TheObject)
                fObjectInList <- TRUE
  
    fObjectInList
}
