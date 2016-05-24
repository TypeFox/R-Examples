.infoPrint <- function(details, limit=1, ...,  prefix='.'){
  if(details>=limit){
    cat(paste(rep(prefix, limit), collapse=''), "") 
    cat(...)  
  }
}

.infoPrint2 <- function(details, limit=1, fmt, ...,  prefix='.'){
  if (details>=limit)
    cat(paste(paste(rep(prefix, limit), collapse=''),  sprintf(fmt, ...), collapse=' '))
}

.colstr <- function(x, collapse=" ")
  paste(x, collapse=collapse)
