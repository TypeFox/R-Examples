to.uniform <- function(ref, val=NA){
   the.ecdf <-  ecdf(ref)
   if(length(val)==1 && is.na(val)){
      ans <- the.ecdf(ref)
   } else {
      ans <- the.ecdf(val)
   } 
   return(ans)
}

