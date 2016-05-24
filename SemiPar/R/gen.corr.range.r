########## S-function: gen.corr.range ##########

# For correcting the range of a mean estimate
# in a generalized model.

# Last changed: 27 FEB 2003

gen.corr.range <- function(mean.est,family)
{
   if (family=="binomial")
   {
      if (min(mean.est)<=0)
      {
         wt <- (max(mean.est) - 0.05)/(max(mean.est) - min(mean.est))
         mean.est <- (1-wt)*max(mean.est) + wt*mean.est
      }

      if (max(mean.est)>=1)
      {
         wt <- (0.95 - min(mean.est))/(max(mean.est) - min(mean.est))
         mean.est <- (1-wt)*min(mean.est) + wt*mean.est
      }
   }

   if (family=="poisson")
   {
      if (min(mean.est)<=0) # Transform to [1,max(mean.est)]
      {
         wt <- (max(mean.est) - 1)/(max(mean.est) - min(mean.est))
         mean.est <- (1-wt)*max(mean.est) + wt*mean.est
      }
   }

   return(mean.est)
}

######### End of S-function gen.corr.range ##########

