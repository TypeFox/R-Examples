penalty.constraint.LNLN <-
function(dataset,parameter,smoothingpar=10^5) {
   # this function should only be called if TY>1
   m <- ncol(dataset)
   TY <- length(parameter)/(m+1)
   mu <- parameter[TY:((m+1)*TY-1)]
   sigma <- parameter[(m+1)*TY]
      
   pen <- 0
   # build a matrix such that the g.th column contains the mu values for gene g
   mu <- matrix(mu,byrow=T,ncol=m)
   # for all genes 1,....m and all types 1,...,TY-1:
   for (g in 1:m) {
      max.datapoint <- max(dataset[,g])
      for (i in 1:(TY-1)) {
         # mode of lognormal i
         mode.ig <- exp(mu[i,g]-sigma^2)
         if (mode.ig<max.datapoint) {
            # sequence of values where to compare the densities
            x <- seq(mode.ig,max.datapoint,(max.datapoint-mode.ig)/500)
            # add mode of population i+1 to this sequence if it is on the right of the mode of population i
            mode.igplus1 <- exp(mu[i+1,g]-sigma^2)
            if (mode.igplus1>mode.ig) {
               x <- c(x,mode.igplus1)
               x <- x[order(x)]
            }             
            # density of lognormal i at these values
            d.LN1 <- d.sum.of.lognormals(y=x,mu.vector=mu[i,g],sigma.vector=sigma)
            # density of lognormal i+1 at these values
            d.LN2 <- d.sum.of.lognormals(y=x,mu.vector=mu[i+1,g],sigma.vector=sigma)
            # penalize if d.LN2>d.LN1
            pen <- pen + sum(pmax(0,d.LN2-d.LN1)^2)
            # now check whether population i+1 is peaked; in that case, the function d.sum.of.lognormals 
            # most probably smoothed the density too much
            if (sum(d.LN2)==0) {
               if (mode.igplus1>mode.ig) {
                  pen <- pen + max(0,dlnorm(mode.igplus1,mu[i+1,g],sigma)-dlnorm(mode.igplus1,mu[i,g],sigma))^2
               }
            }            
         }
      }
   }
   return(smoothingpar*pen)
}
