toSymbols <- function(x, ...)
  UseMethod("toSymbols")

toSymbols.integer <- function(x, ...){
  if(!is.integer(x))
    stop("'x' must be an integer")
  
  symb <- elements[match(x, elements[, "num"]), "symb"]
  symb <- as.character(symb)
  
  if(any(is.na(x)))
    warning("NAs introduced by coercion")
  
  return(symb)
}

toSymbols.numeric <- function(x, ...){
  if(!is.numeric(x))
    stop("'x' must be an numeric")
  
  if(any(round(x) != x))
    stop("'x' must be a whole number")
  
  symb <- toSymbols(as.integer(x))
  
  return(symb)
}

toSymbols.character <- function(x, nletters = 3, ...)
{
  x <- sub(' +$','',sub('^ +', '', x))
  x <- gsub("[0-9]","",x)

  l1 <- substr(x, 1,1)
  ln <- substr(x, 2, nletters)
  
  x <- paste0(toupper(l1), tolower(ln))

  x <- elements[match(x, elements[,"symb"]),"symb"]
  x <- as.character(x)

  if(any(is.na(x)))
    warning("NAs introduced by coercion")
  
  return(x)
}