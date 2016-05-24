countCDFxt <- function(mu, sig2, nt, Nc, Ne, tq=nt, tmax=50, Nboot=500, plot=TRUE)
{
   SEmu <- sqrt(sig2/tq)       # calculate standard error of mu
   # initialize array to store lower and upper bootstrap confidence limits
   # for extinction probabilities
   Glo <- matrix(1,tmax,1)
   Gup <- matrix(0,tmax,1)
   # Calculates confidence interval for mu using equation 3.12 on p. 70
   CI_mu <- c(mu - qt(0.975, nt-1)*SEmu, mu + qt(0.975, nt-1)*SEmu)
   # Calculates confidence interval for sigma^2
   CI_sig2 <- c((nt-1)*sig2/qchisq(.975,nt-1),(nt-1)*sig2/qchisq(.025,nt-1))
   # Calculate the extinction time CDF for the best parameter estimates
   Gbest <- extCDF(mu,sig2,Nc,Ne,tmax)
   Gbest <- matrix(Gbest,tmax,1)
   # Calculate bootstrap confidence limits for extinction probabilites
   for (i in 1:Nboot)
   {
      # Generate random variates of mu and sig2 within their confidence intervals
      ##  CHANGED murnd <- mu + SEmu * rnorm(1,mu,sqrt(sig2))
      murnd<-Inf
      while (murnd < CI_mu[1] | murnd > CI_mu[2])
      {
         # random number from normal distribution:
         ##  CHANGED murnd <- mu + SEmu * rnorm(1,mu,sqrt(sig2))
         murnd <- mu + SEmu * rnorm(1)
      } 
      # CHANGED sig2rnd <- sig2*rchisq(1,nt-1)/(nt-1)
      sig2rnd<-Inf
      while (sig2rnd < CI_sig2[1] | sig2rnd > CI_sig2[2])
      {
         # random number from chi-square distribution
         sig2rnd <- sig2*rchisq(1,nt-1)/(nt-1)
      } 
      # Calculate extintion probabilities given murnd and sig2rnd
      G <- extCDF(murnd,sig2rnd,Nc, Ne,tmax)
      #Store extreme values in each bootstrap interation
      # (only if value is more extreme than previous bootstrap iterations)
      for (x in 1:tmax)
      {
         if (G[x] < Glo[x]){ Glo[x] <- G[x] }
         if (G[x] > Gup[x]){ Gup[x] <- G[x] }
      }
   } # end bootstrap
   # plot CDF
   if(plot)
   {
     plot(Gbest, log="y", type='l', pch=16, col="blue",
       ylim=c(min(Gbest[Gbest!=0]), max(Gup)),          
       main="Extinction CDF",
       xlab="Years into the future",
       ylab="Cumulative probability of quasi-extinction")   
     lines(Glo, col="red", lty=2)
     lines(Gup, col="red", lty=2)
    }
    data.frame(Gbest,Glo,Gup)
} 

