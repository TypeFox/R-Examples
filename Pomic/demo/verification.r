oask <- devAskNewPage(dev.interactive(orNone = TRUE))
oldpar <- par(mfrow=c(3,3))
#Verification with real Metropolis algorithm for model with likelihood function (without POMDEV calculations)
###########################
#Prepare the field pattern#
###########################
nUpd<-10000
nsamples<-500
set.seed(1313)#set the same random seed as in "demoPOMIC.r" to obtain the same "field distribution pattern"
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
#Preparation of the storage container
vec0<-matrix(NA,nrow=nUpd,ncol=2) #2 columns: 'a' and deviance values
for(n in 1:nUpd)
{
  #propose candidate values for parameters
  candidate_a<-rnorm(1,actual_a,sigma_kernel)
  #calculate prior's 
  prior_actual_a<-dnorm(actual_a,mu_a,sigma_prior_a) 
  prior_candidate_a<-dnorm(candidate_a,mu_a,sigma_prior_a)
  #calculate likelihoods
  likelihood_candidate<-dnorm(field,candidate_a,fixed_sd)
  likelihood_actual<-dnorm(field,actual_a,fixed_sd)
  #compute the probability to accept the new parameter set 
  hp<-(prior_candidate_a/prior_actual_a)*prod((likelihood_candidate/likelihood_actual))
  if(runif(1) < hp)
  {
    #candidate is accepted:
    actual_a<-candidate_a
  }
  #store the value of the chain and save the deviance measure
  deviance<- -2*sum(log(likelihood_candidate))
  vec0[n,]<-c(actual_a,deviance)
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
vec<-matrix(NA,nrow=nUpd,ncol=2) 
for(n in 1:nUpd)
{
  candidate_a<-rnorm(1,actual_a,sigma_kernel)
  prior_actual_a<-dnorm(actual_a,mu_a,sigma_prior_a) 
  prior_candidate_a<-dnorm(candidate_a,mu_a,sigma_prior_a)
  likelihood_candidate<-dnorm(field,candidate_a,fixed_sd)
  likelihood_actual<-dnorm(field,actual_a,fixed_sd)
  hp<-(prior_candidate_a/prior_actual_a)*prod((likelihood_candidate/likelihood_actual))
  if(runif(1) < hp)
  {
    actual_a<-candidate_a
  }
  deviance<- -2*sum(log(likelihood_candidate))
  vec[n,]<-c(actual_a,deviance)
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
#Preparation of the storage container
vec2<-matrix(NA,nrow=nUpd,ncol=3) #3 columns: 'a', 'b' and deviance value
for(n in 1:nUpd)
{
  #propose candidate values for parameters
  candidate_a<-rnorm(1,actual_a,sigma_kernels[1])
  lb_b<-ifelse(actual_b-sigma_kernels[2]>0,actual_b-sigma_kernels[2],0)
  ub_b<-ifelse(actual_b+sigma_kernels[2]<range_prior_b,actual_b+sigma_kernels[2],range_prior_b)
  candidate_b<-runif(1,lb_b,ub_b)
  #check if candidate should be accepted
  prior_actual_a<-dnorm(actual_a,mu_a,sigma_prior_a) 
  prior_candidate_a<-dnorm(candidate_a,mu_a,sigma_prior_a)
  likelihood_candidate<-dnorm(field,candidate_a,candidate_b)
  likelihood_actual<-dnorm(field,actual_a,actual_b)
  likratio<-prod((likelihood_candidate/likelihood_actual))
  hp<-(prior_candidate_a/prior_actual_a)*ifelse(is.finite(likratio),likratio,0)   #here we do not need the prior of b since it is uniform
  if(runif(1) < hp)
  {
      actual_a<-candidate_a
      actual_b<-candidate_b
  }
  deviance<- -2*sum(log(likelihood_candidate))
  vec2[n,]<-c(actual_a,actual_b,deviance)
}
plot(as.ts(vec2[,1]),main="History of parameter a, model2")
plot(as.ts(vec2[,2]),main="History of parameter b, model2")
plot(density(vec2[,1]),main="Posterior distribution of a, model2")
plot(density(vec2[,2]),main="Posterior distribution of b, model2")

#We observe that the three model's chains are stable and can thus be used as estimate of posterior
#We can now calculate the two expected deviance measure values for each model and the DIC:
#model0
MeanDEV0<-mean(vec0[,2])
EstimateDEV0<- -2*sum(log(dnorm(field,mean(vec0[,1]),1)))
pD0<-MeanDEV0 - EstimateDEV0
DIC0<-2*MeanDEV0 - EstimateDEV0
#model1
MeanDEV1<-mean(vec[,2])
EstimateDEV1<- -2*sum(log(dnorm(field,mean(vec[,1]),fixed_sd)))
pD1<-MeanDEV1 - EstimateDEV1
DIC1<-2*MeanDEV1 - EstimateDEV1
#model2
MeanDEV2<-mean(vec2[,3])
EstimateDEV2<- -2*sum(log(dnorm(field,mean(vec2[,1]),mean(vec2[,2]))))  
pD2<-MeanDEV2 - EstimateDEV2
DIC2<-2*MeanDEV2 - EstimateDEV2

#preparing a final table showing the results of calculations of DIC
final<-data.frame(model0=c(mean(vec[,1]),NA,MeanDEV0,EstimateDEV0,pD0,DIC0),model1=c(mean(vec[,1]),NA,MeanDEV1,EstimateDEV1,pD1,DIC1),model2=c(mean(vec2[,1]),mean(vec2[,2]),MeanDEV2,EstimateDEV2,pD2,DIC2))
rownames(final)<-c("Estimate a","Estimate b","MeanDEV","EstimateDEV","pD","DIC")
print(final)

par(oldpar)
devAskNewPage(oask)
