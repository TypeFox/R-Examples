#count the number of occurrences of a substring within a string
options(encoding="utf-8")
str_count <- function(x, pattern, sep=""){
  unlist(lapply(
    strsplit(x, sep),
    function(z) na.omit(length(grep(pattern, z)))
  ))
}
#check list current, if it doesn't exist in the list then add this element
add_unique<-function(list,value) {
  if(!is.null(value))
  {
    ifelse (value %in% list, return(list) , return(c(list,value)) ) 
  }
}
#show all file un the result
xfile <- function(sep=":")
{
  tryCatch(
  {
  conf = fromJSON(paste(.libPaths()[1], "x.ent/www/config/ini.json", sep='/'))
  data <- readLines(conf$result$file)
  lst <- c();
  #read line to line
  for(i in 1:length(data))
  {
    if(nchar(data[i]) > 0)
    {
      v <- unlist(strsplit(data[i], sep))[1]
      lst <- add_unique(lst,v)
    }
  }
  return(lst)
  },
  error=function(cond) {
    message("There are problems in paths, please use command xconfig() for verifying your parameters!")
    return(NA)
  },
  warning=function(cond) {
    message("There are problems in paths, please use command xconfig() for verifying your parameters!")
    return(NULL) 
  },
  finally={
    rm(list=ls())
  })
}

