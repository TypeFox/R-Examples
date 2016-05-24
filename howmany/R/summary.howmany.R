"summary.howmany" <-
function(object,...)
  {
    m <- object$m
    m0 <- length(object$pvalues)
    lower <- max(object$lowerbound)
    maxabs <- which.max(object$lowerbound)
    switch(min(3,lower+1),
           {
         
             cat("\n Multiple testing of ",m," hypotheses. \n \n At confidence level ",1-object$alpha, ", no evidence for false null hypotheses was found.\n \n ",sep="")
           },
           {
             cat("\n Multiple testing of ",m," hypotheses. \n \n At confidence level ",1-object$alpha, ", there is at least ", lower," correct rejection, \n among the first ", maxabs  ," rejections. \n \n ",sep="")
             
           },
           {
             cat("\n Multiple testing of ",m," hypotheses. \n \n At confidence level ",1-object$alpha, ", there are at least ", lower," correct rejections \n (all among the first ", maxabs  ," rejections). \n \n ",sep="")
             proportion <- object$lowerbound/(1:m0)
             maxprop <- max((1:m0)[proportion==max(proportion)])
             if( maxprop>0 &  (object$lowerbound[maxprop]/maxprop)>1.2* (lower/maxabs)  )
               {
                 if(maxprop>1)
                   {
                     cat("Furthermore, at confidence level ",1-object$alpha, ", there are at least ",object$lowerbound[maxprop] ," correct rejections among the first ", maxprop," rejections. \n \n ",sep="")
                   }else{
                     cat("At confidence level ",1-object$alpha, ", the first rejection is furthermore correct. \n \n ",sep="")
                   } 
               }
           }
           )
    cat("\n")
  }

