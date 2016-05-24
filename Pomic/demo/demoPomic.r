############################################################
##                                                        ##
##   Demo of POMDEV, metropolis sampling                  ##
##       and POMIC calculations, as found in Appendix B   ##
##       of Piou et al. (2009)			                      ##
##                                                        ##
############################################################
require(Pomic)
oask <- devAskNewPage(dev.interactive(orNone = TRUE))
oldpar <- par(mfrow=c(3,3))
#Metropolis type of algorithm for stochatic model with pseudo-likelihood estimated with POMDEV
###########################
#Prepare the field pattern#
###########################
nUpd<-10000
nsamples<-500
set.seed(1313)#set a random seed (to memorize for verification)
field<-rnorm(nsamples,1,1) 

#########################
#       Model 0         #
#########################
#model 0 is a model generating a normal distribution but we need to find a parameter determining the mean ('a'). The standard deviation is fixed to 1 in this model and should thus become exactly as the original data
fixed_sd<-1
#we assume a loose normal prior on a ~ N(0,1E3)
sigma_prior_a<-1000
actual_a<-mu_a<-0 
##we assume a gaussian kernel for 'a' with size following the rule of thumb (assume 'a' should be found between -10 and 10)
sigma_kernel<-20/sqrt(nUpd)
##Preparation of a start value for an evolutive POMDEV value that will evolve as the conservative pseudo-likelihood measure to be compared to when proposing new parameter values.
modelini<-rnorm(nsamples,mu_a,fixed_sd)
POMDEVevo<-pomdev(field,modelini)
#Preparation of the storage container
vec0<-matrix(NA,nrow=nUpd,ncol=2) #2 columns: 'a' and POMDEV scores
for(n in 1:nUpd)
{
  #propose candidate values for parameters
  candidate_a<-rnorm(1,actual_a,sigma_kernel)
  #create model results
  model<-rnorm(nsamples,candidate_a,fixed_sd)
  #calculate POMDEV_prime
  POMDEV_prime<-pomdev(field,model)
  #check if candidate should be accepted
  prior_actual_a<-dnorm(actual_a,mu_a,sigma_prior_a)  #probability to have actual 'a' value from prior
  prior_candidate_a<-dnorm(candidate_a,mu_a,sigma_prior_a)  #probability to have proposed 'a' value from prior
  likelihood_ratio<-exp((POMDEVevo-POMDEV_prime)/2) #calculate pseudo-likelihood ratio from the POMDEV scores of the mean of past accepted values (POMDEVevo) and the new one
  hp<-(prior_candidate_a/prior_actual_a)* likelihood_ratio #compute the probability to accept the new parameter set
  if(runif(1) < hp)
  {
    #candidate is accepted:
    actual_a<-candidate_a
    #change the conservative POMDEVevo value so that strong likelihood increase do not block the MCMC to temporal (and eventually local) minimum due to the stochastic aspect of the model creating the data
    POMDEVevo<-(POMDEV_prime+(POMDEVevo*(n-1)))/n
  }
  #store the value of the chain and save the deviance measure
  vec0[n,]<-c(actual_a,POMDEV_prime)
}
#plot the history and parameter posterior distribution (if the chain is stable)
plot(as.ts(vec0[,1]),main="History of parameter a, model0")
plot(density(vec0[,1]),main="Posterior distribution of parameter a, model0")

#########################
#       Model 1         #
#########################
#model 1 is a model identical to model 0 but the standard deviation is fixed to 2 in this model. This model will be a bit harder to adjust because of this higher variance.
fixed_sd<-2
sigma_prior_a<-1000
actual_a<-mu_a<-0 
sigma_kernel<-20/sqrt(nUpd)
modelini<-rnorm(nsamples,mu_a,fixed_sd)
POMDEVevo<-pomdev(field,modelini)
vec<-matrix(NA,nrow=nUpd,ncol=2)
for(n in 1:nUpd)
{
  candidate_a<-rnorm(1,actual_a,sigma_kernel)
  model<-rnorm(nsamples,candidate_a,fixed_sd)
  POMDEV_prime<-pomdev(field,model)
  prior_actual_a<-dnorm(actual_a,mu_a,sigma_prior_a) 
  prior_candidate_a<-dnorm(candidate_a,mu_a,sigma_prior_a)
  hp<-(prior_candidate_a/prior_actual_a)*exp((POMDEVevo-POMDEV_prime)/2)
  if(runif(1) < hp)
  {
    actual_a<-candidate_a
    POMDEVevo<-(POMDEV_prime+(POMDEVevo*(n-1)))/n
  }
  vec[n,]<-c(actual_a,POMDEV_prime)
}
plot(as.ts(vec[,1]),main="History of parameter a, model1")
plot(density(vec[,1]),main="Posterior distribution of parameter a, model1")

#########################
#       Model 2         #
#########################
#model 2 is a model generating also a normal distribution but we need to find the mean (a) and standard deviation (b) as parameter for this model
#we assume a uniform prior on 'b' ]0;20] and a loose normal prior on 'a' ~ N(0,1E3)
sigma_prior_a<-1000
range_prior_b<-20   
actual_a<-mu_a<-0 
actual_b<-mu_b<-2
#we assume a gaussian kernel for 'a' and a uniform for 'b'
#we set the size of the kernel a bit wider for 'a' than for the previous models because of the influence of searching also for 'b'
sigma_kernels<-c(40/sqrt(nUpd),range_prior_b/(sqrt(nUpd)))
#Preparation of a start value for the evolutive POMDEV value
modelini<-rnorm(nsamples,mu_a,mu_b)
POMDEVevo<-pomdev(field,modelini)
#Preparation of the storage container
vec2<-matrix(NA,nrow=nUpd,ncol=3) #3 columns: 'a', 'b' and POMDEV scores
for(n in 1:nUpd)
{
  #propose candidate values for parameters
  candidate_a<-rnorm(1,actual_a,sigma_kernels[1])
  lb_b<-ifelse(actual_b-sigma_kernels[2]>0,actual_b-sigma_kernels[2],0)
  ub_b<-ifelse(actual_b+sigma_kernels[2]<range_prior_b,actual_b+sigma_kernels[2],range_prior_b)
  candidate_b<-runif(1,lb_b,ub_b)
  #create model results
  model<-rnorm(nsamples,candidate_a,candidate_b)
  #calculate POMDEV_prime
  POMDEV_prime<-pomdev(field,model)
  #check if candidate should be accepted
  prior_actual_a<-dnorm(actual_a,mu_a,sigma_prior_a) 
  prior_candidate_a<-dnorm(candidate_a,mu_a,sigma_prior_a)
  hp<-(prior_candidate_a/prior_actual_a)*exp((POMDEVevo-POMDEV_prime)/2)    #here we do not need the prior of b since it is uniform
  if(runif(1) < hp)
  {
    actual_a<-candidate_a
    actual_b<-candidate_b
    POMDEVevo<-(POMDEV_prime+(POMDEVevo*(n-1)))/n 
  }
  vec2[n,]<-c(actual_a,actual_b,POMDEV_prime)
}
plot(as.ts(vec2[,1]),main="History of parameter a, model2")
plot(as.ts(vec2[,2]),main="History of parameter b, model2")
plot(density(vec2[,1]),main="Posterior distribution of a, model2")
plot(density(vec2[,2]),main="Posterior distribution of b, model2")

#We observe that the three model's chains are stable and can thus be used as estimate of posterior
#We can now calculate the two expected POMDEV values for each model and the POMIC :
#model0
MeanPOMDEV0<-mean(vec0[,2])    #Mean deviance value
EstimateModel0<-rnorm(nsamples,mean(vec0[,1]),1)  #run model0 for the estimate of 'a' from posterior
EstimatePOMDEV0<-pomdev(field,EstimateModel0) #Deviance value with estimated parameter(s) (here the mean of 'a')
pD0<-MeanPOMDEV0 - EstimatePOMDEV0            #Effective number of parameters
POMIC0<-2*MeanPOMDEV0 - EstimatePOMDEV0       #POMIC calculation
#model1
MeanPOMDEV1<-mean(vec[,2])
EstimateModel1<-rnorm(nsamples,mean(vec[,1]),fixed_sd)  #run model1 for the estimate of 'a' from posterior
EstimatePOMDEV1<-pomdev(field,EstimateModel1)
pD1<-MeanPOMDEV1 - EstimatePOMDEV1
POMIC1<-2*MeanPOMDEV1 - EstimatePOMDEV1
#model2
MeanPOMDEV2<-mean(vec2[,3])
EstimateModel2<-rnorm(nsamples,mean(vec2[,1]),mean(vec2[,2]))  #run model2 for the estimates of 'a' & 'b' from posteriors
EstimatePOMDEV2<-pomdev(field,EstimateModel2)
pD2<-MeanPOMDEV2 - EstimatePOMDEV2
POMIC2<-2*MeanPOMDEV2 - EstimatePOMDEV2

#plotting the field data distribution and the model results with parameters as estimates from the posterior
hist(field,freq=F)
lines(density(EstimateModel0),col=2)
lines(density(EstimateModel1),col=3)
lines(density(EstimateModel2),col=4)
legend(x="topleft",legend=c("model0","model1","model2"),col=c(2:4),lty=c(1,1,1))

#preparing a final table showing the results of calculations of POMIC
final<-data.frame(model0=c(mean(vec0[,1]),NA,MeanPOMDEV0,EstimatePOMDEV0,pD0,POMIC0),model1=c(mean(vec[,1]),NA,MeanPOMDEV1,EstimatePOMDEV1,pD1,POMIC1),model2=c(mean(vec2[,1]),mean(vec2[,2]),MeanPOMDEV2,EstimatePOMDEV2,pD2,POMIC2))
rownames(final)<-c("Estimate a","Estimate b","MeanPOMDEV","EstimatePOMDEV","pD","POMIC")
print(final)

par(oldpar)
devAskNewPage(oask)
