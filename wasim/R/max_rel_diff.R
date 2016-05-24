   max_rel_diff <- function(x,y){
      return(max(diff(x) / diff(y), na.rm=TRUE))
   }

