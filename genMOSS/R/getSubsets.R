getSubsets <-
function (index, tools) {

  if (index == 1) 
    return (1) 
  else { 
    subsetList <- currentWorkLoad <- na.omit(tools$downLinks[index,])
    for (j in 1:length(currentWorkLoad))
      subsetList <- c(subsetList, getSubsets(currentWorkLoad[j], tools))
    return (unique(subsetList))
  }

}
