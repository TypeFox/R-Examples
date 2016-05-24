dfRemoveCases <- function(Data,Cases)
#Returns dataframe with Cases (identified by rowname) removed.  generates warning if cases not found
{ 
  if (is.numeric(Cases))  
  {
    d = Data[!is.element(as.numeric(rownames(Data)),Cases), ]
    FoundIndex =  is.element(Cases,as.numeric(rownames(Data)))
    if (sum(FoundIndex) < length(Cases))
    {      
      warning(c('The following cases were not found in dataframe: ', paste(setdiff(Cases,as.numeric(rownames(Data))), collapse= ', ')))
    } 
  }  else
  {
      d = Data[!is.element(rownames(Data),Cases), ]
      FoundIndex =  is.element(Cases,rownames(Data))
      if (sum(FoundIndex) < length(Cases))
      {       
        warning(c('The following cases were not found in dataframe: ', paste(setdiff(Cases,rownames(Data)), collapse= ', ')))
      } 
  }  
  return(d)
}