bart = function(
   x.train, y.train, x.test=matrix(0.0,0,0),
   sigest=NA, sigdf=3, sigquant=.90, 
   k=2.0,
   power=2.0, base=.95,
   binaryOffset=0,
   ntree=200,
   ndpost=1000, nskip=100,
   printevery=100, keepevery=1, keeptrainfits=TRUE,
   usequants=FALSE, numcut=100, printcutoffs=0,
   verbose=TRUE
)
{
   binary=FALSE

   if(is.factor(y.train)) {
      if(length(levels(y.train)) != 2) stop("y.train is a factor with number of levels != 2")
      binary = TRUE
      y.train = as.numeric(y.train)-1
   } else {
      if((length(unique(y.train)) == 2) & (max(y.train) == 1) & (min(y.train) == 0)) {
         cat('NOTE: assumming numeric response is binary\n')
         binary = TRUE
      }
   }

   if(is.vector(x.train) | is.factor(x.train)) x.train = data.frame(x=x.train)
   if(is.vector(x.test) | is.factor(x.test)) x.test = data.frame(x=x.test)

   if(is.data.frame(x.train)) {
      if(nrow(x.test)) {
         if(!is.data.frame(x.test)) stop('x.train is a data frame so x.test must be also')
	 xtemp = rbind(x.train,x.test)
	 xmtemp = makeind(xtemp)
	 x.train = xmtemp[1:nrow(x.train),,drop=FALSE]
	 x.test = xmtemp[nrow(x.train) + 1:nrow(x.test),,drop=FALSE]
      } else {
         x.train = makeind(x.train)
      }
   }


   #check input arguments:
   if((!is.matrix(x.train)) || (typeof(x.train)!="double")) stop("argument x.train must be a double matrix")
   if((!is.matrix(x.test)) || (typeof(x.test)!="double")) stop("argument x.test must be a double matrix")
   if(!binary) {
      if((!is.vector(y.train)) || (typeof(y.train)!="double")) stop("argument y.train must be a double vector")
   }
   if(nrow(x.train) != length(y.train)) stop("number of rows in x.train must equal length of y.train")
   if((nrow(x.test) >0) && (ncol(x.test)!=ncol(x.train))) stop("input x.test must have the same number of columns as x.train")
   if((!is.na(sigest)) && (typeof(sigest)!="double")) stop("input sigest must be double")
   if((!is.na(sigest)) && (sigest<0.0)) stop("input sigest must be positive")
   if((mode(sigdf)!="numeric") || (sigdf<0)) stop("input sigdf must be a positive number")
   if((mode(printevery)!="numeric") || (printevery<0)) stop("input printevery must be a positive number")
   if((mode(keepevery)!="numeric") || (keepevery<0)) stop("input keepevery must be a positive number")
   if((mode(sigquant)!="numeric") || (sigquant<0)) stop("input sigquant must be a positive number")
   if((mode(ntree)!="numeric") || (ntree<0)) stop("input ntree must be a positive number")
   if((mode(ndpost)!="numeric") || (ndpost<0)) stop("input ndpost must be a positive number")
   if((mode(nskip)!="numeric") || (nskip<0)) stop("input nskip must be a positive number")
   if((mode(k)!="numeric") || (k<0)) stop("input k must be a positive number")
   if(mode(numcut)!="numeric") stop("input numcut must be a numeric vector")
   if(length(numcut)==1) numcut = rep(numcut,ncol(x.train))
   if(length(numcut) != ncol(x.train)) stop("length of numcut must equal number of columns of x.train")
   numcut = as.integer(numcut)
   if(min(numcut)<1) stop("numcut must be >= 1")
   if(typeof(usequants)  != "logical") stop("input usequants must a logical variable")
   if(typeof(keeptrainfits)  != "logical") stop("input keeptrainfits must a logical variable")
   if(typeof(verbose)  != "logical") stop("input verbose must a logical variable")
   if(mode(printcutoffs)  != "numeric") stop("input printcutoffs must be numeric")
   printcutoffs = as.integer(printcutoffs)
   if(printcutoffs <0) stop("input printcutoffs must be >=0")
   if(power <= 0) stop("power must be positive")
   if(base <= 0) stop("base must be positive")

   rgy = range(y.train)
   y = -.5 + (y.train-rgy[1])/(rgy[2]-rgy[1])

   # if sigest=NA, fit a lm to training data to get the value of sigest...
   # sigest is on the scale of the transformed y, so we do the lm after the scaling above...
   if(!binary) {
      if (is.na(sigest)) {
         templm = lm(y~x.train) 
         sigest = summary(templm)$sigma
      } else {
         sigest = sigest/(rgy[2]-rgy[1]) #put input sigma estimate on transformed scale
      }
   } else {
      sigest=1
   }

   ncskip = floor(nskip/keepevery)
   ncpost = floor(ndpost/keepevery)
   nctot = ncskip + ncpost
   totnd = keepevery*nctot

   if(!binary) {
      binaryOffset = -1000.0
   }

   cres = .C('mbart',as.integer(nrow(x.train)), as.integer(ncol(x.train)), as.integer(nrow(x.test)),
                   as.double(x.train), as.double(y),
                   as.double(x.test),
                   as.double(sigest),   as.integer(sigdf), as.double(sigquant),
                   as.double(k),
		   as.double(power), as.double(base),
		   as.double(binaryOffset),
		   as.integer(ntree),      as.integer(totnd),
                   as.integer(printevery), as.integer(keepevery),  as.integer(keeptrainfits),
                   as.integer(numcut), as.integer(usequants), as.integer(printcutoffs),
		   as.integer(verbose),
                   sdraw=double(nctot), 
                   trdraw=double(nrow(x.train)*nctot),
                   tedraw=double(nrow(x.test)*nctot),
                   vcdraw=integer(ncol(x.train)*nctot))
	
   # now read in the results...
   if(!binary) {
   sigma = cres$sdraw*(rgy[2]-rgy[1])
   first.sigma = sigma[1:ncskip] # we often want the sigma draws 
   sigma = sigma[ncskip+(1:ncpost)]

   # put sigest on the original y scale for output purposes
   sigest = sigest*(rgy[2]-rgy[1]) 
   }

   yhat.train = yhat.test = yhat.train.mean = yhat.test.mean = NULL
   varcount = NULL

   if (keeptrainfits) {
      yhat.train = matrix(cres$trdraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]
      if(!binary) {
         yhat.train = (rgy[2]-rgy[1])*(yhat.train+.5) + rgy[1]
         yhat.train.mean = apply(yhat.train,2,mean)
      }
   }
   if (nrow(x.test)) {
      yhat.test = matrix(cres$tedraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]
      if(!binary) {
         yhat.test = (rgy[2]-rgy[1])*(yhat.test+.5) + rgy[1]
         yhat.test.mean = apply(yhat.test,2,mean)
      }
   }

   if(binary) {
      if(keeptrainfits) yhat.train = yhat.train + binaryOffset
      if(nrow(x.test)) yhat.test = yhat.test + binaryOffset
   }

   varcount = matrix(cres$vcdraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]

   if(binary) {
   retval = list(
      call=match.call(),
      yhat.train=yhat.train,
      yhat.test=yhat.test,
      varcount=varcount,
      binaryOffset = binaryOffset
   )
   } else {
   retval = list(
      call=match.call(),
      first.sigma=first.sigma,
      sigma=sigma,
      sigest=sigest,
      yhat.train=yhat.train,
      yhat.train.mean=yhat.train.mean,
      yhat.test=yhat.test,
      yhat.test.mean=yhat.test.mean,
      varcount=varcount,
      y = y.train
   )
   }
   class(retval) = 'bart'
   return(invisible(retval))
}
