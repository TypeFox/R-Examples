rand.test <-
function(set1, set2, sims=1000, crit=.95, graph=TRUE, seed=2) {
   set1 <- data.frame(set1)
   set2 <- data.frame(set2)

   # This part gets the data ready to be used in the randomization test
   samp.distr=c()     #Create a sampling distribution vector for the average absolute r
   samp.distsig=c()   #Create a sampling distribution vector for the number statistically significant
   complete = complete.cases(cbind(set1,set2)) # Combine the data sets and keep only complete cases
   set1.set = subset(set1, subset=complete) #Store the "complete" data sets
   set2.set = subset(set2, subset=complete)
   n = nrow(set1.set) #Find the sample size
   critT = qt(.025, n-2, lower.tail=FALSE) # Find the critical t (assumes alpha = .05) for each test
   critr = sqrt( critT^2 / (critT^2 + n - 2) ) # Find the critical r value
   AbsRObs =  mean(abs(cor(set1.set, set2.set))) #Find the Avg. Absolute R Obs
   SigObs = sum(abs(cor(set1.set, set2.set)) >= critr) #Find the number significant observed

  # This part starts the randomization
  if(seed!=F) {set.seed(seed)}
    for (i in 1:sims) {
     rand.order = sample(n, n, replace=FALSE)   #Generate a sample of random orders
     cor.mat = cor(set1.set[rand.order,],set2.set)  #Get the simulated correlation matrix
     samp.distr[i] = mean(abs(cor.mat))  #Store the absolute average simulated r's in samp.distr
     samp.distsig[i] = sum(abs(cor.mat) >= critr) #Store the number significant in samp.distsig
    }

   # This part computes the statistical properties of the two sampling distributions
   SimMeanR = mean(samp.distr)  #Compute the mean of the sampling distribution
   SimSDr = sd(samp.distr)    #And the SD
   Crit95r = quantile(samp.distr,crit)     #And the critical value (default 95th percentile)
   pr = sum(samp.distr >= AbsRObs) / sims #Find the probability of the observed value
   pr.me <- sqrt(pr * (1-pr) / sims) * qnorm(.9995) 
   SimMeanSig = mean(samp.distsig) #Compute the mean
   SimSDsig = sd(samp.distsig)     # SD
   Crit95Sig = quantile(samp.distsig,crit)  # Critical value (default 95th percentile)
   pSig = sum(samp.distsig >= SigObs) / sims  # Compute a probability value
   pSig.me <- sqrt(pSig * (1-pSig) / sims) * qnorm(.9995)
   if(pr + pr.me > 1.00 | pr - pr.me < .00) {warning("Confidence intervals for p-values may be inaccurate. Try a larger number of sims.")}
   if(pSig + pSig.me > 1.00 | pSig - pSig.me < .00) {warning("Confidence intervals for p-values may be inaccurate. Try a larger number of sims.")}
   if(pr == .00 | pSig == .00) {warning("When p=.00, confidence interval for p not valid.")}
   

   #Clean up and print the results
   out.AbsR = round(rbind(n, AbsRObs, SimMeanR, SimSDr, pr, pr + pr.me, pr - pr.me, Crit95r),4)
   colnames(out.AbsR) = c("Average Absolute r")
   rownames(out.AbsR) = c("N", "Observed", "Exp. By Chance", "Standard Error", "p", "99.9% Upperbound p", "99.9% Lowerbound p", "95th %")
   out.Sig = round(rbind(n, SigObs, SimMeanSig, SimSDsig, pSig, pSig + pSig.me, pSig - pSig.me, Crit95Sig),4)
   colnames(out.Sig) = c("Number Significant")
   rownames(out.Sig) = c("N", "Observed", "Exp. By Chance", "Standard Error", "p", "99.9% Upperbound p", "99.9% Lowerbound p", "95th %")
   results <- list("AbsR"=out.AbsR, "Sig"=out.Sig)

   # This part creates histogram graphics of the sampling distributions
   if (graph == TRUE) {
    old.par = par(mfrow=c(2,1))  # Sets the PAR command two produce two vertical histograms
    hist(samp.distr, freq=TRUE, col="cyan",      #Create a histogram of the sampling distribution
      main="Approximate Sampling Distribution \n For Average Absolute r",
      xlab = "Average Absolute r", ylab="Frequency",
      xlim= range(min(samp.distr)-.01,AbsRObs+.01) )
    abline (v=(Crit95r), col="red")           #Plot the critical value as a line
    points(AbsRObs,0, col="red", pch=19)     #Plot the observed value point
   hist(samp.distsig, freq=TRUE, col="cyan",      #Create a histogram of the sampling distribution
      main="Approximate Sampling Distribution \n For Number Significant",
      xlab = "Number Statistically Significant", ylab="Frequency",
      xlim= range(min(samp.distsig)-1,(SigObs+1)))
    abline (v=(Crit95Sig), col="red")           #Plot the critical value as a line
    points(SigObs,0, col="red", pch=19)     #Plot the observed value point
   }
  return(results)
}
