#Calculate the span
span <- function(x){
      toRet <- if( sum(is.na(x)) == length(x)){
          NA
      } else {
          diff(range(x,na.rm=TRUE))
      }

         if(toRet==-Inf){
            cat("diff_range = -Inf !")
            browser()
         }
      return(toRet)
   }

