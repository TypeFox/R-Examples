##                    author.search                      ##
##      This code is part of the timetree package        ##
## copyright F.-S. Krah 2015 (last update: 2015-04-15)   ##
author.search <- function(author){
  url <- paste("http://www.timetree.org/time_query.php?author_query=",
               author, "&submit=Search", sep="")
  parse <- readHTMLTable(url)
  br <- grep(intToUtf8(194), parse[[1]][,1])  
  f <- vector("list", length(br))
  for(i in seq_along(br)){
    f[[i]] <- rep(i, br[1])
    if(i>1){
      f[[i]] <- rep(i, br[i]-br[i-1])
    }
  }
  f <- unlist(f)
  pubs <- split(parse[[1]], f)
  pubs <- lapply(pubs, na.omit)
  obj <- list(publ.nr = length(pubs), publ=pubs)
  return(obj)
}