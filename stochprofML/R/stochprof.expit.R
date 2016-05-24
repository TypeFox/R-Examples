stochprof.expit <-
function(y,u=0,v=1) {
# The inverse logit function, which maps numbers y from the real line to the interval [u,v]:
#
# expit(y) = u + (v-u) * exp(y)/(1+exp(y)).
#                    
# This function is used here instead of the provided default function
# because there are some numerical instabilities which cause errors in
# our application.             

   large <- (y>=500)
   small <- (y<=(-1000))   
   
   res <- u + (v-u)*exp(y)/(1+exp(y))
   res[large] <- v
   res[small] <- u
   
   return(res)   
}
