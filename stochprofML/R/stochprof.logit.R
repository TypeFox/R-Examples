stochprof.logit <-
function(x,u=0,v=1) {
# The logit function, which maps numbers x from the interval [u,v] to the real line:
#
# logit(x) = log((x-u)/(v-x))     for x in [u,v].
#                    
# This function is used here instead of the provided default function
# because there are some numerical instabilities which cause errors in
# our application.                     
                        
   large <- which(round(x,8)==v)
   small <- which(round(x,8)==u)
   together <- c(large,small)
   
   if ((round(x,8)<u) || (round(x,8)>v)) {
      stop("stochprof.logit: x is outside the allowed parameter range.")
   }   
   
   if (length(together)>0) {
      res <- rep(NA,length(x))
      res[-together] <- log((x[-together]-u)/(v-x[-together]))
      res[large] <- 500
      res[small] <- -1000
   }
   else {
      res <- log((x-u)/(v-x))
   }

   return(res)
}
