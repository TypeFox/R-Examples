remove.dup.rows <- function(dfr) {
  o <- do.call("order",dfr)
#  isdup <- do.call("cbind",lapply(dfr[o,],function(x) eql(x,c(x[-1],NA))))
#  all.dup <- apply(isdup, 1, all)
  all.dup <- do.call("pmin",lapply(dfr[o,],function(x) eql(x,c(x[-1],NA))))
  all.dup[o] <- all.dup 
  dfr[!all.dup,]
}

##  i.e. sort the dataframe, figure out which rows have all values
##  identical to their successor. This gives logical vector, but in the
##  order of the sorted values, so reorder it. Finally select nondups. As
##  a "bonus feature", I think this will also remove any row containing all
##  NA's... 
##  
##  A major stumbling block is that you'll want two NAs to compare equal,
##  hence the eql() function.
##  
##  Actually, I think you can do away with the isdup array and do
##  
##  all.dup <- do.call("pmin",lapply(dfr[o,],function(x)eql(x,c(x[-1],NA))))
##  
##  and there may be further cleanups possible.
##  
##  One dirty trick which is much quicker but not quite as reliable is
##   
##  dfr[!duplicated(do.call("paste",dfr)),]
##  
##  (watch out for character strings with embedded spaces and underflowing
##  differences in numeric data!)
##  
