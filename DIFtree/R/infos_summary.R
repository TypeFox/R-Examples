infos_summary <-
function(splits,
                          items,
                          model,
                          type){
  
  nitems <- length(items)
  output <- data.frame("item"=numeric(nitems),
                       "dif"=character(nitems),
                       "type"=character(nitems),
                       "variables"=character(nitems),
                       "nosplits"=numeric(nitems),
                       stringsAsFactors=FALSE)
  
  for(i in 1:nitems){
    info <- info_summary(splits,items[i],model,type)
    output[i,"item"]            <- info$item
    output[i,"dif"]             <- info$dif
    output[i,"type"]            <- info$type
    output[i,"variables"]       <- info$variables
    output[i,"nosplits"]    <- info$nos
  }
  
  return(output)
}
