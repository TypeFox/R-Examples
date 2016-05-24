TimeToTargetRunSize <- function(U, time, targetU, file_prefix, ci_prob=.95){
#
#  Take the results from the MCMC runs and use it to estimate the
#  time until the target value of the cumulative U is found. For example,
#  what is the time when 20,000 fish have passed through the system.
# 
#  Arguments:
#    U = sim.list from the MCMC run. Usually results$sims.list$U
#    time = time values for strata. Usually results$data$time
#    targetU = target value
#    file_prefix = file prefix for plot
#    ci_prob = what level of credible interval do you want?
# 
#  Returns posterior sample of TimeToTargetRunSize along with summary statistics and quantiles
#  Plots the posterior and adds some basic summary statistics on the plot
#
#  This will plot the posterior distribution, put the summary values
#  of the target distribution on the plot, and return the usual summary
#  statistics (mean, median, and credible intervals, and sample from
#  the posterior).

index <- rep(0,nrow(U))  #array to save times to reach the target value

for(i in 1:nrow(U)){
   cumU <- cumsum(U[i,])
   T <- approx( cumU, 1:ncol(U), xout=targetU, rule=2)  # find when the target value is reached
   index[i] <- T$y
}
index <- index + min(time) - 1 # convert to strata units

mean.index <- mean(index)
med.index  <- median(index)
sd.index   <- sd(index)
probs <- c(seq(0,1,.05),round((1-ci_prob)/2,3),round(1-(1-ci_prob)/2,3))
quant.index <- quantile(index, probs, names=TRUE)  # get the quantiles

# Generate the density plot and give the relevant statistics as well
pdf(file=paste(file_prefix,"-target.pdf",sep=""))
temp<- density(index)

plot(temp,
    main=paste("Posterior distribution of time needed to reach U=",targetU),
    xlab="Stratum")
text(min(index),max(temp$y), 
    label=paste("Mean  : ",         round(mean.index,1), 
                ";       SD: ",     round(sd.index,1),
                ";       ",round(100*ci_prob,0),"% CI: ", 
                round(quantile(index,prob=(1-ci_prob)/2),1), 
                round(quantile(index,prob=1-(1-ci_prob)/2),1)),pos=4)
dev.off()

return( list(targetU=targetU, mean=mean.index, median=med.index, sd=sd.index, quantiles=quant.index, index=index))

}

