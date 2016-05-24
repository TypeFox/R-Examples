##########################################################################################
#######   TLE-FLXmclust - MvtNormal model driver                                ##########
#######                                                                         ##########
#######   Density function according current parameter estimate.                ##########
#######                                                                         ##########
#######   data  - model data, expected to be matrix     ;                       ##########
#######   estim - model estimate.                                               ##########
#######           Current parameters extracted via coef function;               ##########
#######                                                                         ##########
##########################################################################################

### ... is used for user-supplied stuff and to enable TLE to call Density with 4 paramters
FLXmclust.Density =function(data,estim,model,...) {
  mat <- data[[1]]
  ## require(mvtnorm)
  pars    = coefmclust(estim)
  n       = dim(mat)[1]                # number of cases, changed from data to mat
  nc      = length(pars)                # number of components/clusters
                                        #
  lc=ll   = matrix(data = 0, nrow = n, ncol = nc) # loglikelihood by clusters
  prior   = rep(0,nc)

  ## remove old stuff that exptected model.frame
  ## model   = attr(data,"terms")          # model formula
  ## attributes(model)=NULL
  ## model   = as.formula(model) 
  ## changed because if different input parameters...
  ##mat = model.matrix(model,data)        # design matrix
  ## mat     = as.matrix(data) #model.matrix(model,data)                    # design matrix
                                        #
  for (i in 1:nc)
    {       tmp      = try(dmvnorm(mat, pars[[i]]$center, pars[[i]]$cov)) # density
            if(inherits(tmp,"try-error")) return(NULL)
            prior[i] = pars[[i]]$prior
            ll[,i]   = tmp*prior[i]
            lc[,i]  = tmp
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
  lik     = rowSums(ll)
  if(nc>1)lik[na] = NA
  ##
  return(list(ll=ll,lc=lc,c=c,cc=cc,lik=lik))
}
