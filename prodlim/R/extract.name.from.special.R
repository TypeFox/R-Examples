extract.name.from.special <- function(x,pattern="[()]"){
  if (length(x)==1) 
    rev(unlist(strsplit(x,pattern)))[1]
  else
    as.character(sapply(x,extract.name.from.special))
}
