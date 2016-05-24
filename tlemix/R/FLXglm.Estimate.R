#########################################################################################
#######	  TLE-flexmix gaussian, poisson & binomial model driver			#########
#######										#########
#######	  Returns flexmix class estimate.					#########
#######	  data	- model data, expected to be model.frame;			#########
#######	  ind	- data subset;							#########
#######										#########
#######	  ...	- other parameters of current estimate (GLM) procedure.         #########
#######										#########
#######   ntry = maximum number of trials
##########################################################################################

flexmix.Estimate=function(data,ind=NULL,
  nc,class="hard",cluster=NULL,niter=200,minprior=0.1,model=NULL,family=NULL,ntry=9) {
  ##	classification=c("auto", "weighted", "hard", "random")
  ##require(flexmix)                      
                                        
  ## removed getting family as attribute because it is supplied as argument
  ## set default family
  if(is.null(family)) family = "gaussian" # gaussian, poisson and binomial family

  ## model needs to be supplied
  ## if it's not supplied, stop
  if(is.null(model)) {
    stop("No formula supplied.")
  }
                                        
  if(!is.null(ind)) {
    data=data[ind,]                     # substr attribute is not appropriate
    if(!is.null(cluster)) cluster=cluster[ind,]
  } else data=data[,]
  if(!is.null(cluster)) nc=dim(cluster)[2]
  nobs=dim(data)[1]
                                        #
  control	= list(iter=niter,tol=1.e-9,class=class,verb=0,minprior=minprior)

  ## maximum number of trials, now user supplied, default 9
  ## ntry = 9
  itry = 0								# current try number
                                        
  while(itry<=ntry){
    nmix <-	try(flexmix(model,                 
                        data=data.frame(data), 
                        k=nc,                  
                        cluster=cluster, # =NULL	or cluster=d$c,
                        control=control, # control=list(verb=verb,iter=iter,tol=tol,class=class)
                        model=FLXglm(model, family=family)),silent=TRUE)
    if(!inherits(nmix, "try-error")) break
    if(itry>ntry) break                 # no more than ntry tries
    itry<-itry+1
  }
  if(itry>ntry) return(NULL)            # no solution
  stop=attr(nmix,"stop")                # ready
  if(is.null(stop)) stop=F

  if(T&&!stop){ ### 
	if(class!="weighted" && class!="auto")
      {	control$class="weighted"
		control$iter=nmix@iter      # no more iteration than in previous step
                                        #control$iter=3
		tmp_nmix<-nmix

		nk = nmix@k                     # number of components
		c  = nmix@cluster               # cluster
		cc = NULL                       # Problem: it needs to be matrix
		for(i in 1:nk) {q=c; q[q!=i]=0; cc=cbind(cc,q)} # clusters as matrix
		cc[cc!=0]=1

		old.logLik = nmix@logLik
		itry =0                         # current try number
		while(itry<=ntry)
          {	#Sys.sleep(0.01)					# needed for synchronisation of	the other windows
			nmix <-	try(flexmix(model,  # as.matrix(d)~1,
                                data=data.frame(data),
                                k=nk,       #
                                cluster=cc, # =Previous estimate
                                control=control, # control=list(verb=verb,iter=iter,tol=tol,class=class)
                                model=FLXglm(model, family=family)),silent=TRUE)
			if(!inherits(nmix, "try-error")) break
			if(itry>ntry) break         # no more than ntry tries
			itry<-itry+1
          }
		if(itry>ntry) return(NULL)      # no solution
		new.logLik = nmix@logLik
        
        ## if(new.logLik>old.logLik) cat(">")
        ## else if(new.logLik<old.logLik) cat("\<")
        ## else cat("=")
        
      }
  } #..............................................................................
  nc    =nmix@k                   # number of the	components in the estimate
  if(nc==0) return(NULL)          # no solution

  return(nmix)
}
###
### ######################################################################################
### Returns list of coefficents for glm-class models
coefglm=function(nmix,family)
{
	nc = length(nmix@prior)							# number of the	components in the estimate
	if(nc==0) return(NULL)							# no solution
	coef = rep(list(NULL),nc)
	param =	attr(nmix@components[[1]][[1]],"parameters")
	coef = param$coef
	if(family=="gaussian") sigma = param$sigma
	prior =	nmix@prior
	if(nc>1) for( i	in 2:nc)	{
		param=attr(nmix@components[[i]][[1]],"parameters")
		coef = cbind(coef,param$coef)
		if(family=="gaussian") sigma = c(sigma,param$sigma)
	}
	coef = as.matrix(coef)
	colnames(coef) = paste("C",1:nc,sep="")
	names(prior) = colnames(coef)
	if(family=="gaussian") names(sigma) = colnames(coef)

	if(family=="binomial") return(list(coefficients=coef,prior=prior))
	if(family=="poisson")  return(list(coefficients=coef,prior=prior))
	if(family=="gaussian") return(list(coefficients=coef,sigma=sigma,prior=prior))
	return(NULL)
}

