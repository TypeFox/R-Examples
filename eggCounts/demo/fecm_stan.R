##################################################################
## Illustrate faecal egg count modelling
##################################################################
# load libraries
library('eggCounts')
require('coda',quietly = TRUE)
require('rstan',quietly = TRUE)
##################################################
##  1-sample situation
##################################################
# simulate faecal egg counts
set.seed(123)
counts <- simData1s(n = 20,  # number of samples samples
                    mean = 500, # mean 
                    kappa = 1, # overdispersion parameter
                    phi=1,   # prevalence
                    f=50)    # correction factor

# look at simulated counts: matrix with columns 
# (obs = raw counts*correction factor, master= raw counts ,true = true EpG)
head(counts,3)

# run a Stan model
resultsP <- fec_stan(counts[,"obs"],zeroInflation = FALSE)

# get samples and covert them to mcmc-objects, so that the functionality 
# of the R-package coda (summary, plot,...) can be used
# get samples 

result_mcmc<-stan2mcmc(resultsP)
# this is a list with
#  fec  - mean(EpG rate)

# for instance
plot(result_mcmc)

##################################################
##  2-sample situation
##################################################
# load epgs data before and after treatment
data(epgs)
# plot FECs
epgsL <- reshape(epgs, direction="long",varying=list(c("before","after")))
epgsL$time <- factor(epgsL$time, levels=1:2, labels=c("untreated","after treatment"))

if (require('lattice'))
  xyplot(before/50 ~ time, group=id, data=epgsL, type=c("p","l"), col=1,
         xlab="", ylab="Faecal egg counts")

# run a paired zero-inflation model
result2 <- fecr_stan(epgs$before, epgs$after, preCF=10, paired = TRUE, zeroInflation = TRUE)

# get samples 
result_mcmc2<-stan2mcmc(result2)
# this is a list with
#  fecr=reduction in means),
#  mean(EpG rate in untreated animals) and
#  mean(EpG rate in treated animals).

# plot samples
plot(result_mcmc2[,1:3])
