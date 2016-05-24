##########################################################################################
#######   TLE-FLXmclust - MvtNormal model driver                                 #########
#######                                                                         #########
#######   Returns flexmix class estimate.                                       #########
#######   data  - model data, expected to be model.frame;                       #########
#######   ind   - data subset;                                                  #########
#######                                                                         #########
#######   ...   - other parameters of current estimate (GLM) procedure.         #########
#######                                                                         #########
#########################################################################################

FLXmclust.Estimate=function(data,ind=NULL,
  nc,class="hard",cluster=NULL,niter=200,minprior=0.1,model=NULL,ntry=9,...) {
  ##       classification=c("auto", "weighted", "hard", "random")
  ## require(flexmix)
  ##
  if(is.null(model)) {                  # extract model formula
    ## model = attr(data,"terms")
    ## attributes(model)=NULL
    ## model = as.formula(model)
    stop("Model formula must be supplied.")
  }
  ##
  ## create model matrix, but data must be a matrix, not a data.frame
  ## TODO: ensure that this is the case!
  ## mat=data
  ## mat = as.matrix(data)                 # user-supplied data frame
  var_name <- names(data) # extract variable name
  mat <- data[[1]] # get data from list
  if(!is.null(ind)) {
    mat=mat[ind,]                       # substr attribute is not appropriate
    if(!is.null(cluster)) cluster=cluster[ind,]
  }

  else mat=mat[,]

  if(!is.null(cluster)) nc=dim(cluster)[2]
  nobs=dim(mat)[1]

  if(is.null(nobs)) return(NULL)

  control	= list(iter=niter,tol=1.e-9,class=class,verb=0,minprior=minprior)

  ## build data list again
  data <- list(mat)
  names(data) <- var_name

  ## moved this to function call, should be user supplied
  ##  ntry = 9
  ## maximal number of trials
  itry = 0                              # current try number

  while(itry<=ntry) {
    nmix <- try(flexmix(model, #mat~1,                                      # as.matrix(d)~1,
                        data=data,
                        k=nc,            #
                        cluster=cluster, # =NULL or cluster=d$c,
                        control=control, # control=list(verb=verb,iter=iter,tol=tol,class=class)
                        model=FLXmclust(diagonal=FALSE)),silent=TRUE)
    if(!inherits(nmix, "try-error")) break
    if(itry>ntry) break                 # no more than ntry tries
    itry<-itry+1
  }
	
			
			
	
  if(itry>ntry) return(NULL)            # no solution
  stop=attr(nmix,"stop")                # ready
  if(is.null(stop)) stop=F

  ##................................................................................
  if(T&&!stop) { # the last iterations must be carried out with posterior weights!
	if(class!="weighted" && class!="auto") {
      		control$class="weighted"
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
		while(itry<=ntry) {
          ##nmix <-	try(flexmix(mat~1,				# as.matrix(d)~1,
          nmix <- try(flexmix(model,
                              data,
                              k=nk,       #
                              cluster=cc, # =Previous estimate
                              control=control, # control=list(verb=verb,iter=iter,tol=tol,class=class)
                              model=FLXmclust(diagonal=FALSE)),silent=TRUE)
        if(!inherits(nmix, "try-error")) break
        if(itry>ntry) break             # no more than ntry tries
        itry<-itry+1
      }
    }
    if(itry>ntry) return(NULL)          # no solution
    new.logLik = nmix@logLik
    ## if(new.logLik>old.logLik) cat(">")
    ## else if(new.logLik<old.logLik) cat("\<")
    ## else cat("=")
  }
                                        #..............................................................................
  nc    =nmix@k                 # number of the	components in the estimate
  if(nc==0) return(NULL)        # no solution

  return(nmix)
}

###
#########################################################################################
### # FLXmclust coef function

coefmclust <- function(nmix) {
  nc = length(nmix@prior)          # number of the components in the estimate
  if(nc==0)
    return(NULL)                        # no solution
  coef    = rep(list(NULL),nc)
  names(coef) = paste("C",1:nc,sep="")
  for( i in 1:nc) {
    par     = parameters(nmix,component=i,model=1,simplify=FALSE)
    size    = nmix@size[i]
    names(size) = NULL
    coef[[i]] = list(center=par$center,cov=par$cov,prior=nmix@prior[i],size=size)
  }
  return(coef)
}
