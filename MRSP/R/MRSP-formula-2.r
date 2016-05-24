################################################################################
#####    Functions for handling formulas and model matrices in MRSP.       #####
################################################################################
#####    Author: Wolfgang Pößnecker                                        #####
#####    Last modified: 10.11.2014, 23:52                                  #####
################################################################################

modelmatrix.MRSP <- cmpfun(function(formula, data, modelname){
 ## create the model frame.
 formula <- as.Formula(formula)
 myframe <- model.frame(formula, data)

 ## a little case distinction
 mncase <- modelname == "Multinomial Logit Model"
 cumcase <- modelname == "Cumulative Logit Model"
 seqcase <- modelname == "Sequential Logit Model"

 #### create the response object 'y'
 y <- as.matrix(model.part(formula, myframe, lhs=1))
 if(!(all(y %in% c(0,1)) | all(y %in% c(F,T))))
    stop("The response must be supplied as a 0-1-vector or a logical vector!")
    
 if(ncol(y) > 1){          # this case actually shouldnt appear since we supply data in 'long' format.
  K <- ncol(y)
  nobs <- nrow(y)          # number of individual observations
  ymat <- matrix(0, nrow=nobs, ncol=K)
  ymat[cbind(seq_len(nobs), max.col(y))] <- 1
 }else{
  K <- nrow(y)/sum(y)              # this works since y must be filled with either 0-1 or False-True.
  nobs <- sum(y)
  if(is.logical(y)){
   ymat <- matrix(ifelse(y==T,1,0), nrow=nobs, ncol=K, byrow=T)
  }else  if(all(y %in% c(0,1))){
   ymat <- matrix(y, nrow=nobs, ncol=K, byrow=T)
  }
 }
 if((cumcase | seqcase) && all(rowSums(ymat)==1)) ymat <- ymat[,-K]      # in ordinal regression, we use a response matrix with K-1 columns

 ###### now the covariates part. due to case distinctions, things become a bit messy... ######
 len.form <- length(attr(formula, "rhs"))
 if(mncase && len.form > 3)
    stop("You cannot use more than three different combinations of predictor and effect types in multinomial logit models!")
 if(cumcase && len.form > 1)
    stop("You cannot use anything but global predictors with global effects in cumulative logit models!")
 #if(seqcase && len.form > 4)
 #   stop("you cannot use more than four different combinations of predictor and effect types in sequential logit models!")
 if(seqcase && len.form > 2)
    stop("You can currently only use global predictors, with either global or category-specific coefficients, in sequential logit models!",'\n',
    "  Support for category-specific predictors will be added in a future release of MRSP.",'\n')

 #### create the matrix of global/individual-specific covariates, called 'x'
 x <- model.matrix(formula, myframe, rhs=1)
 ## MRSP currently cant handle cases in which there are no global predictors:
 if(ncol(x) < 1)
   stop("MRSP currently cannot handle models without at least one global predictor. (For example, the column of ones that yields the class-specific intercepts.)")

 if(nrow(x) == nobs*K) x <- x[seq(1,nobs*K, by=K),]            # this ensures that we use only the actual x-observations and no duplicates if the data is supplied in long format
 which.intercept.x <- which(apply(x, 2, function(u) diff(range(u)) < .Machine$double.eps ^ 0.5))
 has.intercept.x <- !(length(which.intercept.x) == 0)
 X <- x
 x.vars <- attr(terms(formula, rhs=1), "term.labels")

 ## if we are using a sequential logit model, there might be two parts of the formula belonging to x:
 ## one for those global covariates with category-/alternative-specific effects and one with global effects.
 if(seqcase && len.form > 1){
  x2 <- model.matrix(formula, myframe, rhs=2)
  if(nrow(x2) == nobs*K) x2 <- x2[seq(1,nobs*K, by=K),]
  ## handle potentially double intercepts
  which.intercept.x2 <- which(apply(x2, 2, function(u) diff(range(u)) < .Machine$double.eps ^ 0.5))
  has.intercept.x2 <- !(length(which.intercept.x2) == 0)
  if(has.intercept.x2){
   x2 <- x2[,-which.intercept.x2]
  }
  X <- cbind(x, x2)          # we need to keep x and x2 for the creation of penindex and grpindex
  x2.vars <- attr(terms(formula, rhs=2), "term.labels")
 }
 
 #### create the object called 'Z' which contains the category-/alternative-specific covariates
 ## first the z-covariates with global/individual-specific effects:
 if((len.form > 1 && mncase) | (len.form > 2 && seqcase)){
  if(mncase){
   z <- model.matrix(formula, myframe, rhs=2)
   z.vars <- attr(terms(formula, rhs=2), "term.labels")
  }else if(seqcase){
   z <- model.matrix(formula, myframe, rhs=3)
   z.vars <- attr(terms(formula, rhs=3), "term.labels")
  }
  ## since intercepts are always global and thus part of x, intercepts in z must be removed
  which.intercept.z <- which(apply(z, 2, function(u) diff(range(u)) < .Machine$double.eps ^ 0.5))
  has.intercept.z <- !(length(which.intercept.z) == 0)
  if(has.intercept.z){
   z <- z[,-which.intercept.z]
  }
  ## create a data object in the form that is required by MRSP.fit:
  Z <- list(); length(Z) <- K
  for(i in seq_len(K)){
   Z[[i]] <- z[seq(i,(nobs*K-(K-i)),by=K), ]
  }
 }
 
 ## now the z-covariates with category-/alternative-specific effects:
 if((len.form > 2 && mncase) | (len.form > 3 && seqcase)){
  z2 <- model.matrix(formula, myframe, rhs=len.form)
  z2.vars <- attr(terms(formula, rhs=len.form), "term.labels")
  ## since intercepts are always global and thus part of x, intercepts in z must be removed
  which.intercept.z2 <- which(apply(z2, 2, function(u) diff(range(u)) < .Machine$double.eps ^ 0.5))
  has.intercept.z2 <- !(length(which.intercept.z2) == 0)
  if(has.intercept.z2){
   z2 <- z2[,-which.intercept.z2]
  }
  ## create a data object in the form that is required by MRSP.fit:
  for(i in seq_len(K)){
   Z[[i]] <- cbind(Z[[i]], z2[seq(i,(nobs*K-(K-i)),by=K), ])
  }
 }

 #### prepare the output
 out <- list(formula=formula, y=ymat, X=X, x=x, nobs=nobs, len.form=len.form, x.vars=x.vars) ##/## which.intercept=which.intercept.x, has.intercept=has.intercept.x, )
 if(exists("x2",inherits=F)) out$x2 <- x2
 if(exists("x2.vars",inherits=F)) out$x2.vars <- x2.vars
 if(exists("Z",inherits=F)) out$Z <- Z
 if(exists("z",inherits=F)) out$z <- z
 if(exists("z.vars",inherits=F)) out$z.vars <- z.vars
 if(exists("z2",inherits=F)) out$z2 <- z2
 if(exists("z2.vars",inherits=F)) out$z2.vars <- z2.vars
 
 return(out)
})




















