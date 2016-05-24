### #############################################################################################
### #############################################################################################
### ####	  TLE.Flexmix binomial,poisson & gaussian model driver			       ##############
### ####										                                   ##############
### ####	  Density function according current parameter estimate.		       ##############
### ####										                                   ############## 
### ####	  data	- model data, expected to be model.frame with family attribute;##############
### ####	  pars	- model estimate.						                       ##############
### ####		  Current parameters extracted via coef function;		           ##############
### ####                                                                           ##############
### #############################################################################################

flexmix.Density	=function(data,pars,model,family) {
  ##
  ## family <- attr(data,"family") ## family is now user supplied..
  ## this should be written in a nicer way - TODO for the next version)
  est	= coefglm(pars,family) # changed to coefglm # get current parameter estimate
  par	= est$coefficients
  sigma	= est$sigma
  n	= dim(data)[1]
  ## nc	= pars@k
  nc	= length(est$prior)             # number of components
  ##
  ##    family	= attr(data,"family")						# model family
  if(is.null(family)) family="guassian"
  ##
  ## model	= attr(data,"terms")						# model formula
  ##	  attributes(model)=NULL
  ## model	= as.formula(model) # now supplied as argument
  mat	= model.matrix(model,data)      # design matrix
                                        #
  ## this is a work-around to get the response.
  ## is there a better way to do this?
  yresp   = model.response(model.frame(model,data))        # responce variable(s)
  linp    = mat%*%par                   # linear predictor
  if(family=="poisson") lambda = exp(linp) # poisson model parameter
  if(family=="binomial") {prob = exp(linp); prob = prob/(1+prob) # binomial model probability
                          success = yresp[,1]   # number of successes
                          size = rowSums(yresp) # size of cell
                        }
                                        #
                                        #.......................................................................................
                                        #
  ll	= matrix(data =	0, nrow	= n, ncol = nc)
                                        #
  for (i in 1:nc)	{
    if(family=="gaussian") tmp = try(dnorm(yresp,linp[,i],sigma[i]))
    if(family=="poisson")  tmp = try(dpois(yresp,lambda[,i]))
    if(family=="binomial") tmp = try(dbinom(success,size,prob[,i]))
    if(inherits(tmp,"try-error")) return(NULL)
    ll[,i]	= tmp
  }
  if(nc>1)na  = apply(ll,1,function(x)any(is.na(x)))
  ##
  ##       Cluster affiliations and Loglikelihoods by cases:
  ##
  c	= apply(ll,1,function(x)sort(x,index.return=TRUE,decreasing=TRUE)$ix[1]) # cluster
  cc  	= NULL
  for(i in 1:nc) {q=c; q[q!=i]=0; cc=cbind(cc,q)} # clusters as matrix
  cc[cc!=0]=1
  ##
  lik	= apply(ll,1,function(x)sum(est$prior*x)) # Loglikelihood
  if(nc>1)lik[na] = NA
  ##
  return(list(ll=ll,c=c,cc=cc,lik=lik))
}
