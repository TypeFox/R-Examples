neffects <-
function(x, term)
{
  n <- 0L
  for(i in 1L:length(x)) {
    ne <- names(x[[i]]$effects)
    if(is.null(ne) || !is.character(term))
      ne <- 1L:length(x[[i]]$effects)
    for(j in 1L:length(term)) {
      if(is.character(term[j])) {
        tmp <- splitme(term[j])
        tmp <- resplit(tmp[tmp != " "])
        take <- NULL
        for(jj in 1:length(ne)) {
          ## if(!is.na(pmatch(tmp, ne[jj])))
          if(length(grep(tmp, ne[jj], fixed = TRUE)))
            take <- c(take, jj)
        }
      } else take <- match(term[j], ne)
      if(length(take) > 0L && length(x[[i]]$effects) > 0L && !is.na(take))
        n <- n + length(take)
    }
  }

  return(n)
}

