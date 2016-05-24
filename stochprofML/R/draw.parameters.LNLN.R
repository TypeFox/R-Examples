draw.parameters.LNLN <-
function(par.range,m) {
# Draws a random parameter from the specified range. The first column of par.range
# contains the lower bounds for all parameters; the second columns the upper bounds.
#
# The so-drawn parameter is modified such that it fulfils the required restrictions:
# - The sum p_1 + ... + p_(T-1) of the probabilities for the first T-1 types 
#   is smaller than or equal to one.
# - The mu-values for the first gene are decreasing.
#
# m is the number of genes analyzed.
#
# There is no parameter "fix.mu". If mu is fixed, the respective lower and upper bounds 
# are identical so that the drawn parameter will be the one to which mu is fixed.
      
   # number of types
   TY <- nrow(par.range)/(m+1)
   
   # if there are infinite bounds, replace them by a finite number that can be 
   # handled by runif
   par.range[par.range==Inf] <- 10^7
   par.range[par.range==-Inf] <- -10^7   
   
   # draw a random parameter vector
   draw <- runif(nrow(par.range),par.range[,1],par.range[,2])

   ###############################################
   # sort the mu values and p values             #
   # such that the mus for gene 1 are decreasing #
   ###############################################
   
   # first, compute the full p vector (the current one has only TY-1 components)
   if (TY>1) {
      p <- draw[1:(TY-1)]
      if (sum(p)>1) {
         p <- c(p/sum(p),0)
      }
      else {
         p <- c(p,1-sum(p))
      }
   }
   else {
      p <- 1
   }
   
   # extract mu for gene 1 and its order
   mu <- draw[TY:(length(draw)-1)]   
   mu <- matrix(mu,byrow=T,ncol=m)
   mu.gene1 <- mu[,1]
   ord <- order(mu.gene1,decreasing=T)
   
   # sort p
   p <- p[ord]
   
   # sort mu
   mu <- mu[ord,]
   
   # plug in back into "draw"   
   # p  
   if (TY>1) { 
      draw[1:(TY-1)] <- p[-length(p)]
   }
   # mu
   draw[TY:(length(draw)-1)] <- c(t(mu))
   
   return(draw)
}
