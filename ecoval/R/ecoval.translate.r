ecoval.translate <- function(word,dictionary)
{
  w <- as.character(dictionary[word])
  if ( is.na(w) ) return(word)
  return(w)
}


ecoval.dict <- function(language,dictionaries=NA)
{
  dict <- NA
  
  if ( !(is.data.frame(dictionaries) | is.matrix(dictionaries)) ) dictionaries <- ecoval::ecoval.dictionaries.default
  
  if ( ! is.na(language) )
  {
    ind <- match(language,colnames(dictionaries))
    if ( ! is.na(ind) )
    {
      dict <- as.character(dictionaries[,ind])
      names(dict) <- dictionaries[,1]
    }
  }
  return(dict)
}
