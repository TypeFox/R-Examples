varRename <- function(Data, From, To)
#2010-10-11: released, JJC
{
  if (length(From) != length(To)) stop("Length of 'From' and 'To' vectors must match")
  
  for (i in 1:length(From) ) {
    if(is.element(From[i], names(Data))){   
      
        if (!is.element(To[i],names(Data))) names(Data)[names(Data)==From[i]] = To[i] 
        else stop('Variable name already exists in dataframe')   
        
    } else warning (sprintf('Variable name does not exist: %s', From[i]))
  }

  return(Data)
}

