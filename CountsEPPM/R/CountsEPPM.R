CountsEPPM <-
function(formula,data,model.type='mean and variance',model='general',
                  offset=NULL,initial=NULL,ltvalue=NA,utvalue=NA,
                  optimization.method='optim',control=NULL,
                  scale.factor.model='no',fixed.b=NA) {

# Checking formula and data

   if (is.data.frame(data)==TRUE) { 
      cat('\n','Dependent variable is a vector of single counts.','\n')
      FBoth     <- Formula(formula) 
      mfBoth    <- model.frame(FBoth,data=data)
      resp.var  <- model.part(FBoth,data=mfBoth,lhs=1)
      nvar      <- length(data)
      nobs      <- length(data[[1]])
      if (nvar==1) { covariates  <- NULL
                   } else { covariates <- data } 
      list.counts  <- lapply(1:nobs, function(i) 
                            c(rep(0,resp.var[[1]][i]),1) )
      mean.obs     <- resp.var[[1]]
      variance.obs <- rep(0,nobs)
      scalef.obs   <- rep(0,nobs)

                                   } else { 

      cat('\n','Dependent variable is a list of frequency distributions of counts','\n')

      nvar      <- length(data)
      nobs      <- length(data[[1]])
      if ((nvar==1) & (is.list(data[[1]])==TRUE)) { list.counts <- data[[1]] 
                                                    covariates  <- NULL
          } else { for ( i in 1:nvar ) { 
            if (is.list(data[[i]])==TRUE)  { list.counts <- data[[i]] 
                            if (i==1) { covariates   <- data.frame(c(rep(NA,nobs)))
                                             } else { 
                            covariates   <- data.frame(covariates,c(rep(NA,nobs))) }
                          } else { if (i==1) { covariates   <- data.frame(data[[1]])
                                             } else { 
                            covariates   <- data.frame(covariates,data[[i]]) }
                                 } } # end for loop
                                                 } # end if nvar==1 & is.list

      if (is.null(covariates)==FALSE) { names(covariates) <- names(data) }
      mean.obs     <- rep(0,nobs)
      variance.obs <- rep(0,nobs)
      scalef.obs   <- rep(0,nobs)
      for ( i in 1:nobs ) { 
         count  <- list.counts[[i]]
         ncount <- sum(count)  
         nmax1  <- length(count)
         nmax   <- nmax1 - 1
         cnum   <- 0:nmax
         if (ncount==1) { mean.obs[i] <- t(cnum)%*%count
                        } else {
            mean.obs[i]     <- t(cnum)%*%count/ncount 
            variance.obs[i] <- (t(cnum*cnum)%*%count - 
                            ncount*mean.obs[i]*mean.obs[i]) / (ncount - 1)
            if (mean.obs[i]>0) {
               scalef.obs[i]   <- variance.obs[i]/mean.obs[i] }
# Checking for NA for variance.obs and scalef.obs introduced Version 1.01
            if (is.na(variance.obs[i])==TRUE) { variance.obs[i] <- 0 }
            if (is.na(scalef.obs[i])==TRUE)   { scalef.obs[i] <- 0 }
                              } } # end of for loop

                            }  # end of if is.data.frame
   wkdata <- data.frame(mean.obs,variance.obs,scalef.obs)
# Checking for no covariates introduced Version 1.01
   if (is.null(covariates)==FALSE) { 
      wkdata <- data.frame(wkdata,covariates) }
# Checking arguments of function

# Setting error indicator
   vinerror <- rep(0,8)

# Checking for correct model.types 
   if ((model.type!='mean only') & (model.type!='mean and variance')) {
      covariates.matrix.mean     <- NULL
      covariates.matrix.variance <- NULL
      vinerror[1] <- 1 }

# if vinerror[1] is 1 there is no point in checking formula

inerror <- sum((vinerror>0))

if (inerror==0) {   

# Checking that formula has 1 lhs and either 1 or 2 rhs
# corresponding to just a formula for mean, and a formula for both mean
# and variance 
   FBoth <- Formula(formula)
   lenFB <- length(FBoth)
   if (model.type=='mean only') {
      FBoth <- update(FBoth,mean.obs ~ .) } # end of mean only
   if (model.type=='mean and variance') {
      if (lenFB[2]==1) { 
         cat('\n','Model for variance or scale factor set to intercept only.','\n')}
         if (scale.factor.model=='yes') { 
            if (lenFB[2]==1) { 
               FBoth <- update(FBoth,mean.obs | scalef.obs ~ . | 1 ) 
                              } else {
               FBoth <- update(FBoth,mean.obs | scalef.obs ~ . | . ) }
                                        } else {
            if (lenFB[2]==1) { 
               FBoth <- update(FBoth,mean.obs | variance.obs ~ . | 1 )
                                        } else {
               FBoth <- update(FBoth,mean.obs | variance.obs ~ . | . ) }
                                        }} # end of if mean and variance 

   wk.nobs          <- length(wkdata[[1]])
   mfBoth <- model.frame(FBoth,data=wkdata)
   lenFB <- length(FBoth)
   covariates.matrix.mean <- model.matrix(FBoth,data=mfBoth,lhs=1,rhs=1)
   resp.mean <- model.part(FBoth,data=mfBoth,lhs=1)
   if (lenFB[1]==1) { covariates.matrix.variance <- 
      matrix(c(rep(1,nrow(covariates.matrix.mean))),ncol=1) } 
   if (lenFB[1]==2) { covariates.matrix.variance <- 
      model.matrix(FBoth,data=mfBoth,lhs=2,rhs=2) } 

                } # end of check of model.type

# Checking for correct models
   if (model.type=='mean only') {
      if ((model!='Poisson') & (model!='negative binomial') & 
          (model!='negative binomial fixed b') & 
          (model!='Faddy distribution') & 
          (model!='Faddy distribution fixed b')) { vinerror[2] <- 1 }}
   if (model.type=='mean and variance') {
      if ((model!='general') & (model!='general fixed b') & 
          (model!='limiting')) { vinerror[2] <- 1 }
# Checking that the formula has two parts 
      if (lenFB[2]==1) {
         cat('\n','variance model set to default of an intercept only','\n') } }

# Checking that scale.factor.model is only yes if model is under- 
# or over-dispersed means & variances  
    if ((scale.factor.model=='yes') & (model.type!='mean and variance')) { 
       vinerror[3] <- 1 }

# Checking that list.counts is a list  
   if (is.list(list.counts)==FALSE) { vinerror[4] <- 1 }

# Checking for offsets i.e. offset.mean, offset.variance are not null
    if (is.null(offset)==TRUE) { 
       offset.mean     <- c(rep(0,nobs)) 
       offset.variance <- c(rep(0,nobs)) 
                                } else {
       if (length(offset)==1) { offset.mean     <- offset[[1]] 
                                offset.variance <- NULL }
       if (length(offset)==2) { offset.mean     <- offset[[1]]
                                offset.variance <- offset[[2]] }
       if (length(offset)>2) { cat('\n','WARNING: offset has more than 2 elements','\n')
                cat('\n','calculations proceed but with the default of 0 for offsets used','\n')
                offset.mean     <- NULL
                offset.variance <- NULL }
       if (is.null(offset.mean)==TRUE) {
             offset.mean <- c(rep(0,nobs)) } else {
             nobswk <- length(offset.mean)
             if (nobswk!=nobs) { offset.mean=c(rep(0,nobs)) 
                cat('\n','WARNING: length of offset.mean not same as number of grouped count vectors','\n')
                cat('\n','calculations proceed but with the default for of 0 offset.mean used','\n')
                               } } # end if offset.mean
       if (is.null(offset.variance)==TRUE) { 
             offset.variance <- c(rep(0,nobs)) } else {
             nobswk <- length(offset.variance)
             if (nobswk!=nobs) { offset.variance=c(rep(0,nobs)) 
                cat('\n','WARNING: length of offset.variance not same as number of grouped count vectors','\n')
                cat('\n','calculations proceed but with the default of 0 for offset.variance used','\n')
                               } } } # end if is.null(offset)
    wkdata <- data.frame(wkdata,offset.mean,offset.variance)

# if any of vinerror are 1 there is no point in checking for initial estimates

inerror <- sum((vinerror>0))

if (inerror==0) {   

# Setting up initial estimates of parameters if not input 
   if (is.null(initial)==TRUE) {
      nmax1   <- length(variance.obs)
      n.var.0 <- sum((variance.obs==0))

# data.frame input or case data input as list, 
# n.var.o is number of variances that are 0
      if (n.var.0==nmax1) {
# switching off warnings from glm usually caused by non integer values for the counts
         options(warn=-1)
# Using glm to obtain initial estimates 
         glm.Results <- glm(formula(FBoth,lhs=1,rhs=1),
                            family=poisson(link="log"),
                            data=wkdata,offset=offset.mean) 
# switching warnings from glm back on
        options(warn=0)                       
                              } else {
# identifying and replacing mean.obs values of 0 with 0.5 because a
# gaussian distribution is used with a log link for the initial estimates
         wkdata$mean.obs <- sapply(1:nobs, function(j) 
              if (wkdata$mean.obs[j]==0) { wkdata$mean.obs[j] <- 0.5
              } else { wkdata$mean.obs[j] <- wkdata$mean.obs[j] } )
# Using glm to obtain initial estimates 
         glm.Results <- glm(formula(FBoth,lhs=1,rhs=1),
                            family=gaussian(link="log"),
                            data=wkdata,offset=offset.mean) 
                                   } # if (n.var.0==nmax1)
      initial.mean        <- coefficients(glm.Results)
      names(initial.mean) <- names(coefficients(glm.Results))

      if (model.type=='mean and variance') {
         if (n.var.0==nmax1) { 
# no observed variance as case not grouped data
# fitting Poisson model for variance to mean as dependent variable
# switching off warnings from glm usually caused by non integer values for the counts
            options(warn=-1)
            glm.Results <- glm(formula(FBoth,lhs=1,rhs=2),
                           family=poisson(link="log"),
                           data=wkdata,offset=offset.variance) 
# switching warnings from glm back on
            options(warn=0)                       
                          } else {      
# identifying and replacing variance.obs values of 0 with 0.5
            if (n.var.0>0) { wkdata$variance.obs <- sapply(1:nobs, function(j) 
                               if (wkdata$variance.obs[j]==0) { wkdata$variance.obs[j] <- 0.5
                               } else { wkdata$variance.obs[j] <- wkdata$variance.obs[j] } )
                       } # end of if (n.var.0>0)          
            glm.Results <- glm(formula(FBoth,lhs=2,rhs=2),
                           family=gaussian(link="log"),
                           data=wkdata,offset=offset.variance) }
         initial.variance        <- coefficients(glm.Results)
         names(initial.variance) <- names(coefficients(glm.Results))
         if (model=='general')  { parameter <- c(initial.mean,initial.variance,0) 
            names(parameter) <- c(names(initial.mean),names(initial.variance),"log(b)") }
         if ((model=='general fixed b') | (model=='limiting')) { 
            parameter <- c(initial.mean,initial.variance) 
            names(parameter) <- c(names(initial.mean),names(initial.variance)) }

# Calculating log likelihood for initial estimates for regression of log(variance) (mean)
# on variance model 
         loglikelihood <- LL.Regression.Counts(parameter,model.type,model,list.counts,
               covariates.matrix.mean,covariates.matrix.variance,
               offset.mean,offset.variance,ltvalue,utvalue,
               scale.factor.model,fixed.b)
#         if (optimization.method=='nlm')   { loglikelihood <- -loglikelihood }

# Calculating log likelihood for initial estimates for log(variance) of 0 for
# the variance model 
         wks <- length(initial.variance)
         if (model=='general')  { wk.parameter <- c(initial.mean,rep(0,wks),0) 
            names(wk.parameter) <- c(names(initial.mean),names(initial.variance),"log(b)") }
         if ((model=='general fixed b') | (model=='limiting')) { 
            wk.parameter <- c(initial.mean,rep(0,wks)) 
            names(wk.parameter) <- c(names(initial.mean),names(initial.variance)) }
         wk.loglikelihood <- LL.Regression.Counts(wk.parameter,model.type,model,list.counts,
               covariates.matrix.mean,covariates.matrix.variance,
               offset.mean,offset.variance,ltvalue,utvalue,
               scale.factor.model,fixed.b)
# Use the parameter estimates from the initial one with largest log likelihood
         if (wk.loglikelihood>loglikelihood) { parameter <- wk.parameter }

                                           } # end of if mean and variance

      if (model.type=='mean only') { 
          if ((model=='Poisson') | (model=="negative binomial fixed b")) {
                                 parameter        <- initial.mean
                                 names(parameter) <- names(initial.mean) }
          if (model=="negative binomial") { parameter <- c(initial.mean,0) 
                           names(parameter) <- c(names(initial.mean),"log(b)") }
          if (model=="Faddy distribution") { parameter <- c(initial.mean,0,0) 
                           names(parameter) <- c(names(initial.mean),"c","log(b)") }
          if (model=="Faddy distribution fixed b") { parameter <- c(initial.mean,0) 
                           names(parameter) <- c(names(initial.mean),"c") }
                                   } # of if model.type mean only
                                   } else {
# Checking length of input initial against model value
      npar.one <- ncol(covariates.matrix.mean)
      numpar <- length(initial) 
      if (model.type=="mean only") { npar <- npar.one 
              loc.c <- npar + 1 
              if ((model=="negative binomial") | 
                  (model=="Faddy distribution fixed b")) { npar <- loc.c }
              if (model=="Faddy distribution") { npar <- npar + 2 } }
      if (model.type=="mean and variance") { 
              npar.two <- ncol(covariates.matrix.variance)
              npar <- npar.one + npar.two 
              if (model=="general")  { npar <- npar + 1 } }
      if (numpar!=npar) { vinerror[5] <- 1 }
      parameter <- initial
      if (is.null(names(initial))==TRUE) { 
         cat('\n','WARNING: initial has no associated names','\n')
         names(parameter) <- 1:numpar } # end if is.null(names)
                           } # end if is.null(initial)

# Checking for improper probability distributions for initial estimates
# Setting up vector vnmax and calculating likelihood for initial estimates
    vnmax  <- sapply(list.counts,length) - 1
    output <- Model.Counts(parameter,model.type,model,covariates.matrix.mean,
              covariates.matrix.variance,offset.mean,offset.variance,scale.factor.model,
              fixed.b,vnmax)
    for ( i in 1:nobs) { probability <- output$probabilities[[i]]
       nmax1 <- vnmax[i] + 1
       wks   <- sum(probability)
       if (is.finite(wks)==TRUE) {
          wks <- round(sum(probability),digits=10)
# rounding error can cause the sum of the probabilities to <0 or >1
# so rounded to 10 decimal places                                        
          if ((wks<0) | (wks>1))  { vinerror[6] <- 1 }
          wks <- sum((round(probability,digits=10)<0) | 
                     (round(probability,digits=10)>1))
          if (wks>0) { vinerror[6] <- 1 }
                                 } else { vinerror[6] <- 1
                                        } # end of if is.finite(wks)
                    } # end of for loop
                } # end of check of initial estimates given first five were okay

# Checking that fixed.b has a positive value if that model is used
    if (((model=='general fixed b') | 
         (model=='negative binomial fixed b distribution fixed b') |
         (model=='Faddy distribution fixed b')) & 
        ((is.na(fixed.b)==TRUE) | (fixed.b<=0))) { vinerror[7] <- 1 }

# Checking that the input initial estimate for c is <1 if a 
# Faddy distribution model is being fitted.
    if ((is.null(initial)==FALSE) & ((model=='Faddy distribution') | 
        (model=='Faddy distribution fixed b'))) {
       if (round(parameter[loc.c],digits=4)>=1) { vinerror[8] <- 1 }}

inerror <- sum((vinerror>0))

##########################################################################
# start of main calculations when inerror==0 
##########################################################################

if (inerror==0) {   

# Checking for truncation i.e. ltvalue and utvalue are not NA
    if (is.na(ltvalue)==FALSE) { cat('\n','distribution truncated below at ',ltvalue) }
    if (is.na(utvalue)==FALSE) { cat('\n','distribution truncated above at ',utvalue,'\n') }

# Setting defaults and checking control parameters for optim and nlm
   if (is.null(control)==TRUE) { 
     if (optimization.method=='optim') { 
        control=list(fnscale=-1,trace=0,maxit=1000,abstol=1e-8,reltol=1e-8,
                     alpha=1.0,beta=0.5,gamma=2.0) 
                                       } # end of if optim
     if (optimization.method=='nlm') { 
        control <- list(fscale=1,print.level=0,stepmax=1,gradtol=1e-12,
                        steptol=1e-14,iterlim=500,ndigit=16)
                                     } # end of if nlm
                              } else { 
# Control parameters related to Nelder-Mead simplex only. Those related
# only to other methods of optim are considered to be unknown and ignored.
     ncontrol <- length(control)
     if (optimization.method=='optim') {
        wk.control <- list(fnscale=-1,trace=0,maxit=1000,abstol=1e-8,reltol=1e-8,
                           alpha=1.0,beta=0.5,gamma=2.0) 
        for ( i in 1:ncontrol) { 
            if ((names(control[i])!='fnscale') & (names(control[i])!='trace') & 
                (names(control[i])!='maxit') & (names(control[i])!='abstol') & 
                (names(control[i])!='reltol') & (names(control[i])!='alpha') & 
                (names(control[i])!='beta') & (names(control[i])!='gamma')) { 
               wkname <- names(control[i])
               cat('\n','WARNING: Argument in control is unknown so ignored.','\n')
                                     } # end of if
            if (names(control[i])=='fnscale') { wk.control$fnscale <- control[[i]] }
            if (names(control[i])=='trace')   { wk.control$trace   <- control[[i]] }
            if (names(control[i])=='maxit')   { wk.control$maxit   <- control[[i]] }
            if (names(control[i])=='abstol')  { wk.control$abstol  <- control[[i]] }
            if (names(control[i])=='reltol')  { wk.control$reltol  <- control[[i]] }
            if (names(control[i])=='alpha')   { wk.control$alpha   <- control[[i]] }
            if (names(control[i])=='beta')    { wk.control$beta    <- control[[i]] }
            if (names(control[i])=='gamma')   { wk.control$gamma   <- control[[i]] }
                                } # end of for loop
        control <- wk.control
                                       } # end of if optim
     if (optimization.method=='nlm') { 
        wk.control <- list(fscale=1,print.level=0,stepmax=1,gradtol=1e-8,
                           steptol=1e-12,iterlim=500,ndigit=12)
        for ( i in 1:ncontrol) { 
            if ((names(control[i])!='fscale') & (names(control[i])!='print.level') & 
                (names(control[i])!='stepmax') & (names(control[i])!='gradtol') &  
                (names(control[i])!='steptol') & (names(control[i])!='iterlim') &
                (names(control[i])!='ndigit')) { 
               wkname <- names(control[i])
               cat('\n','Argument in control is unknown so ignored.','\n')
                                                                      } # end of if
            if (names(control[i])=='fscale')        { wk.control$fscale  <- control[[i]] }
            if (names(control[i])=='print.level')   { wk.control$print.level <- control[[i]] }
            if (names(control[i])=='stepmax')       { wk.control$stepmax <- control[[i]] }
            if (names(control[i])=='gradtol')       { wk.control$gradtol <- control[[i]] }
            if (names(control[i])=='steptol')       { wk.control$steptol <- control[[i]] }
            if (names(control[i])=='iterlim')       { wk.control$iterlim <- control[[i]] }
            if (names(control[i])=='ndigit')        { wk.control$ndigit  <- control[[i]] }
                                } # end of for loop
        control <- wk.control
                                     } # end of if nlm
                                   } # end of is.null(control)

# Fitting model using optim or nlm

        if (optimization.method=='optim') {
           if (length(parameter)==1) {
              Results  <- optim(parameter,fn=LL.Regression.Counts,gr=NULL,model.type,model,list.counts,
                             covariates.matrix.mean,covariates.matrix.variance,
                             offset.mean,offset.variance,ltvalue,utvalue,
                             scale.factor.model,method="Brent",lower=-20,upper=20,hessian=FALSE,
                             control=control) 
                        } else {
              convergence <- list(c("successful"),
                                  c("iteration limit max has been reached"),
                                  c(" "),c(" "),c(" "),c(" "),
                                  c(" "),c(" "),c(" "),c(" "),
                                  c("degeneracy of the Nelder-Mead simplex"))
              Results  <- optim(parameter,fn=LL.Regression.Counts,gr=NULL,model.type,model,list.counts,
                             covariates.matrix.mean,covariates.matrix.variance,
                             offset.mean,offset.variance,ltvalue,utvalue,
                             scale.factor.model,fixed.b,method="Nelder-Mead",hessian=FALSE,
                             control=control)
              wks <- Results$convergence + 1
              cat('\n','optimization method optim:','\n')
              cat(' function calls ',Results$counts[1],'\n')
              cat(' convergence    ',Results$convergence,convergence[[wks]],'\n')
              if (is.null(Results$message)==FALSE) {
                 cat(' message        ',Results$message,'\n') }
                        } # end of if length(parameter)==1
           estimates <- Results$par } # end optimization.method=optim
        if (optimization.method=='nlm') {
           code <- list(c("relative gradient is close to zero,","\n",
                          " current iterate is probably solution","\n"),
                        c("successive iterates within tolerance,","\n",
                          " current iterate is probably solution","\n"),
                        c("last global step failed to locate a point lower than estimate.","\n",
                          " Either estimate is an approximate local minimum of the function,","\n",
                          " the function is too non-linear for this algorithm,","\n",
                          " or steptol is too large","\n"),
                        c("iteration limit max has been reached","\n"),
                        c("Maximum step size stepmax exceeded five consecutive times.","\n",
                          " Either the function is unbounded below, becomes asymptotic to a","\n",
                          " finite value from above in some direction, or stepmax is too small","\n"))
           Results <- nlm(F.Regression.Counts,parameter,model.type,model,list.counts,
                          covariates.matrix.mean,covariates.matrix.variance,
                          offset.mean,offset.variance,ltvalue,utvalue,
                          scale.factor.model,fixed.b,hessian=FALSE,fscale=control$fscale,
                          print.level=control$print.level,stepmax=control$stepmax,
                          gradtol=control$gradtol,steptol=control$steptol,
                          iterlim=control$iterlim,ndigit=control$ndigit)
           cat('\n','optimization method nlm:','\n')
           cat(' iterations  ',Results$iterations,'\n')                     
           cat(' return code ',Results$code,'\n')
           cat(' ',code[[as.numeric(Results$code)]],'\n') 
           estimates <- Results$estimate } # end optimization.method=nlm 
           names(estimates) <- names(parameter) 

          output <- Model.Counts(estimates,model.type,model,covariates.matrix.mean,
                   covariates.matrix.variance,offset.mean,offset.variance,scale.factor.model,
                   fixed.b,vnmax) 

# model.hessian from hessian from numDeriv method Richardson
          model.hessian <- hessian(LL.Regression.Counts,x=estimates,method="Richardson",
                        model.type=model.type,model=model,list.counts=list.counts,
                        covariates.matrix.mean=covariates.matrix.mean,
                        covariates.matrix.variance=covariates.matrix.variance,
                        offset.mean=offset.mean,offset.variance=offset.variance,
                        ltvalue=ltvalue,utvalue=utvalue,
                        scale.factor.model=scale.factor.model,fixed.b=fixed.b)

          npar  <- length(estimates)
# Deleting appropriate row and column of matrix for Faddy distribution if c=1
          if ((model=="Faddy distribution") | 
              (model=="Faddy distribution fixed b")) { nparm1 <- npar - 1
                       nparm1sq <- nparm1*nparm1
                       nparm2   <- nparm1 - 1
                       wk.hessian <- matrix(c(rep(0,nparm1sq)),ncol=nparm1)
            if (model=="Faddy distribution") { wk.c <- round(estimates[nparm1],digits=7)
                          } else { wk.c <- round(estimates[npar],digits=7) }
              if ((model=="Faddy distribution") & (wk.c==1)) { 
                       wk.hessian[1:nparm2,1:nparm2] <- model.hessian[1:nparm2,1:nparm2]
                       wk.hessian[nparm1,] <- c(model.hessian[npar,1:nparm2],model.hessian[npar,npar])
                       wk.hessian[,nparm1] <- c(model.hessian[1:nparm2,npar],model.hessian[npar,npar])
                       model.hessian       <- wk.hessian
                                          } # if Faddy & c=1
              if ((model=="Faddy distribution fixed b") & (wk.c==1)) { 
                       wk.hessian[1:nparm1,1:nparm1] <- model.hessian[1:nparm1,1:nparm1]
                       model.hessian                 <- wk.hessian
                                          } # if Faddy fixed b & c=1
                                                     } # end first if

          deter     <- det(model.hessian)
          if (is.finite(deter)==FALSE) {deter <- 0 }
          condition <- rcond(model.hessian)
          if (is.finite(condition)==FALSE) {condition <- 0 }

# function rcond is from the package Matrix and gives the (reciprocal) condition number
# near 0 is ill-conditioned, near 1 is well conditioned
          if ((abs(deter)<1.e-16) | (abs(condition)<1.e-16)) { se <- rep(NA,npar) 
                        } else { wks <- nrow(model.hessian)
                                 if (wks==1) { se <- sqrt(-1/model.hessian) 
                                  } else { varcov <- - solve(model.hessian) 
                                           se     <- sqrt(diag(varcov)) } } # if end (deter==0)

# Placing appropriate rows of se vector for Faddy distribution if c=1
          if ((model=="Faddy distribution") | 
              (model=="Faddy distribution fixed b")) { 
              if ((model=="Faddy distribution") & (wk.c==1)) { se <- c(se[1:nparm2],NA,se[nparm1])
                                                             } # if Faddy & c=1
              if ((model=="Faddy distribution fixed b") & (wk.c==1)) { se <- c(se[1:nparm1],NA)
                                          } # if Faddy fixed b & c=1
                                                     } # end first if

          loglikelihood <- LL.Regression.Counts(estimates,model.type,model,list.counts,
               covariates.matrix.mean,covariates.matrix.variance,
               offset.mean,offset.variance,ltvalue,utvalue,
               scale.factor.model,fixed.b)
          output.fn <- list(model.type=model.type,model=model,
                       covariates.matrix.mean=covariates.matrix.mean,
                       covariates.matrix.variance=covariates.matrix.variance,
                       offset.mean=offset.mean,offset.variance=offset.variance,
                       ltvalue=ltvalue,utvalue=utvalue,
                       scale.factor.model=scale.factor.model,fixed.b=fixed.b,
                       estses=data.frame(names(estimates),estimates,se),
                       vnmax=vnmax,loglikelihood=loglikelihood,
                       mean.obs=mean.obs,variance.obs=variance.obs)

                 } else {

# Input errors
     vinerror.names <- c('unknown model.type','unknown model for this model.type',
                         'scale factor model should be set to no',
                         'list.counts is not a list','number of parameters error',
                         'improper distribution produced by initial estimates',
                         'value of fixed.b is NA or <=0',
                         'initial c>=1 is not allowed for the Faddy distribution')
     vinerrors <- sapply(1:8, function(i) 
                  if (vinerror[i]==1) { cat('\n',vinerror.names[i],'\n') } )
     output.fn <- list(model.type=model.type,model=model,
                       covariates.matrix.mean=covariates.matrix.mean,
                       covariates.matrix.variance=covariates.matrix.variance,
                       offset.mean=offset.mean,offset.variance=offset.variance,
                       ltvalue=ltvalue,utvalue=utvalue,
                       scale.factor.model=scale.factor.model,fixed.b=fixed.b,
                       estses=data.frame(NA,NA,NA),
                       vnmax=NA,loglikelihood=NA,
                       mean.obs=mean.obs,variance.obs=variance.obs) } # inerror if else
     attr(output.fn, "class") <- c("CountsEPPM")
          return(output.fn)     }
