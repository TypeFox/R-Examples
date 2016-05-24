# 2011-03-10 CJS Utility functions for specifying prior for muTT
 
#****************************************************************************

make.muTT.prior <- function(x){
#
# estimate the muTT prior based on a dirchelete type prior
# x are values representing belief in the travel times.
# For example, x=c(1,4,3,2) represents a system where the
# maximum travel time is 3 strata after release with 
#  1/10=.1 of the animals moving in the stratum of release
#  4/10=.4 of the animals taking 1 stratum to move
#  etc
#
# So if x=c(10,40,30,20), this represent the same movement pattern
# but a strong degree of belief
#
# We convert these into muTT which are the mean(on logit scale) of 
# conditional movement probabilities

   n <- sum(x)
   p <- x/n

   ndelta <-(n-cumsum(c(0,x[-length(x)])))
 
   delta<- x/ndelta
   delta<- delta[-length(delta)] # drop the last term
   sd.delta <- sqrt(delta*(1-delta)/ndelta[-length(ndelta)])

   mean.muTT <- logit(delta)
   sd.muTT <- sd.delta/delta/(1-delta)
 
   list(mean.muTT=mean.muTT, sd.muTT=sd.muTT)
}


#****************************************************************************

visualize.muTT.prior<- function( muTT.prior, npoints=1000, title=NULL){
#
#  Visualize the range of movement rates based on the muTT.prior
#  We will generate normal distribution based on mean and sd and then
#  back transform to actual movement rates which are then plotted
#
#  muTT.prior is a list with 2 element. mean.muTT, and sd.muTT
#   browser()
   p<- matrix(NA, nrow=npoints, ncol=length(muTT.prior$mean.muTT)+1)
   for(i in 1:npoints){
#     generate the movement probability
      muTT <- rnorm(length(muTT.prior$mean.muTT), 
                 mean=muTT.prior$mean.muTT,
                 sd  =muTT.prior$sd.muTT)
#     convert these back to probabilities
      delta <- expit(muTT)
      Theta <- delta[1]
      for(j in 2:length(delta)){
         Theta <- c(Theta, delta[j]*(1-sum(Theta)))
      }
      Theta <- c(Theta,1-sum(Theta))
      p[i,] <- Theta
   }
   #browser()
   colnames(p) <- seq(0,ncol(p)-1,1)
   temp3 <- stack(as.data.frame(p))
   #temp3[1:5,]
   boxplot( values ~ ind, data=temp3,   
      main=if(is.null(title)){'Visualization of prior on movement'} else {title},
      ylab='P(movement)',
      xlab='Time to move')

} # end of function  

