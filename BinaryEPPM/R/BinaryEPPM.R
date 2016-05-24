BinaryEPPM <-
function(formula,data,model.type='p and scale-factor',
                  model='generalized binomial',link='cloglog',
                  offset=NULL,initial=NULL,
                  optimization.method='optim',control=NULL) {

# Checking formula and data

   nvar <- length(data)
   nobs <- length(data[[1]])
   if (nvar==1) { covariates  <- NULL }

   if (is.data.frame(data)==TRUE) { 
      cat('\n','Dependent variable a vector of numerator / denominator.','\n')
      FBoth     <- Formula(formula) 
      mfBoth    <- model.frame(FBoth,data=data)
      resp.var  <- mfBoth[[1]]
      n.var     <- mfBoth[[2]]
      if (nvar>1) { covariates <- data } 
      list.binary  <- lapply(1:nobs, function(i) 
                   (resp.var[i]==0)*(c(1,rep(0,n.var[i]))) +
                   (resp.var[i]==n.var[i])*(c(rep(0,n.var[i]),1)) +
                   ((resp.var[i]>0) & (resp.var[i]<n.var[i]))* 
                   (c(rep(0,resp.var[i]),1,rep(0,(n.var[i]-resp.var[i]))))
                            ) # end of lapply
      mean.obs     <- resp.var
      variance.obs <- rep(0,nobs)
      p.obs <- (resp.var==0)*1.e-14 + 
               ((resp.var>0) & (resp.var<n.var))*resp.var/n.var +
               (resp.var==n.var)*(1-1.e-14) 
      scalef.obs <- rep(0,nobs)
      vnmax      <- n.var

                                   } else { 

      cat('\n','Dependent variable is a list of binomial frequency distributions','\n')
      if ((nvar==1) & (is.list(data[[1]])==TRUE)) { list.binary <- data[[1]] 
          } else { for ( i in 1:nvar ) { 
            if (is.list(data[[i]])==TRUE)  { list.binary <- data[[i]] 
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
      p.obs        <- rep(0,nobs)
      scalef.obs   <- rep(0,nobs)
      vnmax        <- rep(0,nobs)
      for ( i in 1:nobs ) { 
         count        <- list.binary[[i]]
         nmax1        <- length(count)
         nmax         <- nmax1 - 1
         vnmax[i]     <- nmax
         cnum         <- 0:nmax
         ncount       <- sum(count)  
         mean.obs[i]  <- t(cnum)%*%count/ncount 
         p.obs[i]     <- mean.obs[i]/nmax
         if (p.obs[i]==0) { p.obs[i] <- 1.e-10 }
         if (p.obs[i]==1) { p.obs[i] <- 1 - 1.e-10 }
         variance.obs[i] <- (t(cnum*cnum)%*%count - ncount*mean.obs[i]*mean.obs[i]) / (ncount - 1)
         scalef.obs[i]   <- variance.obs[i] / (mean.obs[i]*(1-p.obs[i]))
         if (is.finite(variance.obs[i])==FALSE) { variance.obs[i] <- 0 } 
         if (is.finite(scalef.obs[i])==FALSE) { scalef.obs[i] <- 0 } 
                    } # end of for loop
                            }  # end of if is.data.frame

# Checking for no covariates 
   if (is.null(covariates)==TRUE) {  wkdata <- data.frame(p.obs,scalef.obs) 
               } else { wkdata <- data.frame(p.obs,scalef.obs,covariates) }

# Checking arguments of function

# Setting error indicator
   vinerror <- rep(0,6)

# Checking for correct model.types 
   if ((model.type!='p only') & (model.type!='p and scale-factor')) {
      covariates.matrix.p      <- NULL
      covariates.matrix.scalef <- NULL
      vinerror[1] <- 1 }

# if vinerror[1] is 1 there is no point in checking formula

inerror <- sum((vinerror>0))

if (inerror==0) {   

# Checking that formula has 1 lhs and either 1 or 2 rhs
# corresponding to just a formula for p, and a formula for both p
# and scale-factor 
   FBoth <- Formula(formula)
   lenFB <- length(FBoth)
   if (model.type=='p only') {
      if (lenFB[2]==2) { 
         cat('\n','model.type is p only but rhs of formula has two parts to it')
         cat('\n','2nd part of rhs of formula is ignored','\n') } # end if lenFB[2]
      FBoth <- update(FBoth,p.obs ~ .) } # end of p only
   if (model.type=='p and scale-factor') {
      if (lenFB[2]==1) { 
         cat('\n','Model for scale-factor set to intercept only.','\n')
         FBoth <- update(FBoth,p.obs | scalef.obs ~ . | 1 )
                                  } else {
         FBoth <- update(FBoth,p.obs | scalef.obs ~ . | . ) }
                                         } # end of if p and scale-factor 

   mfBoth  <- model.frame(FBoth,data=wkdata)
   lenFB <- length(FBoth)
   covariates.matrix.p <- model.matrix(FBoth,data=wkdata,lhs=NULL,rhs=1)
   if (lenFB[1]==1) {
      inerror <- 0
      covariates.matrix.scalef <- 
                      matrix(c(rep(1,nrow(covariates.matrix.p))),ncol=1)
                                       } # mean formula only
   if ((lenFB[1]==2) & (lenFB[2]==2)) {
      inerror <- 0
      covariates.matrix.scalef <- model.matrix(FBoth,data=wkdata,lhs=NULL,rhs=2)
                                      } # p and scale-factor formula 

                } # end of check inerror==0

# Checking for correct combinations of model.type and model
   if ((model.type!='p only') & (model.type!='p and scale-factor')) {
      inerror <- 1
      cat('\n','unknown model.type','\n') }
   if (model.type=='p only') {
      if ((model!='binomial') & (model!='generalized binomial') &
          (model!='beta binomial') & (model!='correlated binomial')) {
         inerror <- 1
         cat('\n','unknown model for this model.type','\n') } }
   if (model.type=='p and scale-factor') {
      if ((model!='generalized binomial') & (model!='beta binomial') &
          (model!='correlated binomial')) {
         inerror <- 1
         cat('\n','unknown model for this model.type','\n') }
# Checking that the formula has two parts 
      if ((lenFB[1]==1) & (lenFB[2]==1)) {
         cat('\n','variance model set to default of an intercept only','\n') } }

# Checking that list.binary is a list  
if (is.list(list.binary)==FALSE) { inerror <-1 
   cat('\n','list.binary is not a list','\n') }

   nsuccess <- list.binary
   ntrials  <- list.binary

# Checking for offsets i.e. offset.p, offset.scalef are not null
    if (is.null(offset)==TRUE) { 
       offset.p     <- c(rep(0,nobs)) 
       offset.scalef <- c(rep(0,nobs)) 
                                } else {
       if (length(offset)==1) { offset.p     <- offset[[1]] 
                                offset.scalef <- NULL }
       if (length(offset)==2) { offset.p     <- offset[[1]]
                                offset.scalef <- offset[[2]] }
       if (length(offset)>2) { cat('\n','WARNING: offset has more than 2 elements','\n')
                cat('\n','calculations proceed but with the default of 0 for offsets used','\n')
                offset.p     <- NULL
                offset.scalef <- NULL }
       if (is.null(offset.p)==TRUE) {
             offset.p <- c(rep(0,nobs)) } else {
             nobswk <- length(offset.p)
             if (nobswk!=nobs) { offset.p=c(rep(0,nobs)) 
                cat('\n','WARNING: length of offset.p not same as number of observations','\n')
                cat('\n','calculations proceed but with the default for of 0 offset.p used','\n')
                               } } # end if offset.p
       if (is.null(offset.scalef)==TRUE) { 
             offset.scalef <- c(rep(0,nobs)) } else {
             nobswk <- length(offset.scalef)
             if (nobswk!=nobs) { offset.scalef=c(rep(0,nobs)) 
                cat('\n','WARNING: length of offset.scalef not same as number of observations','\n')
                cat('\n','calculations proceed but with the default of 0 for offset.scalef used','\n')
                               } } } # end if is.null(offset)
    wkdata <- data.frame(wkdata,offset.p,offset.scalef)

# Setting up initial estimates of parameters if not input 
   if (is.null(initial)==TRUE) { 
# data.frame input or case data input as list, 
# Using glm to obtain initial estimates 
# switching off warnings from glm usually caused by non integer values for the counts
      options(warn=-1)
# setting up new formula in standard form for data.frame input
   if (is.data.frame(data)==TRUE) { 
      wkdata <- data.frame(wkdata,resp.var,n.var)
      FBoth_one <- update(FBoth,cbind(resp.var,(n.var-resp.var)) ~ .) 
# changing to standard form of arguments to glm for binomial for data as a data frame
      glm.Results <- glm(formula(FBoth_one,lhs=1,rhs=1),family=binomial(link=link),
                         data=wkdata,offset=offset.p)
                                      } else {
      glm.Results <- glm(formula(FBoth,lhs=1,rhs=1),family=binomial(link=link),
                         data=wkdata,offset=offset.p) } # end if is.data.frame
      initial.p <- coefficients(glm.Results)
# switching warnings from glm back on
      options(warn=0)                       
      if (model.type=='p only') {
          if (model=='binomial') { parameter <- initial.p
                           names(parameter) <- names(initial.p) }
          if (model=='generalized binomial')  { parameter <- c(initial.p,0) 
                           names(parameter) <- c(names(initial.p),"GB parameter") }
          if (model=='beta binomial') { parameter <- c(initial.p,0) 
                           names(parameter) <- c(names(initial.p),"beta binomial theta") }
          if (model=='correlated binomial') { parameter <- c(initial.p,0) 
                           names(parameter) <- c(names(initial.p),"correlated binomial rho") }
          numpar <- length(parameter)
                                   } # of if model.type p only
      if (model.type=='p and scale-factor') {
# switching off warnings from glm usually caused by non integer values for the counts
            options(warn=-1)
# n.var.0 is number of variances that are 0
            n.var.0 <- sum(((scalef.obs==0) | (is.finite(scalef.obs)==FALSE)))
            if (n.var.0==nobs) { 
               glm.Results <- glm(formula(FBoth,lhs=1,rhs=2),family=gaussian(link="log"),
                                  data=wkdata,offset=offset.scalef)              
# the following is done to obtain names
               initial.scalef <- rep(0,length(glm.Results$coefficients))
               names(initial.scalef) <- names(glm.Results$coefficients) 
                               } else {
# observed scale-factors available hence use normal distribution with log link
               if (n.var.0>0) { wk.scalef.obs <- rep(0,nobs)
                                wk.scalef.obs <- sapply(1:nobs, function(i) 
                                   if (scalef.obs[i]==0) { wk.scalef.obs[i] <- 1
                                       } else { wk.scalef.obs[i] <- scalef.obs[i] } ) 
                  FBoth_two <- update(FBoth, . | wk.scalef.obs ~ . | .) 
                  wkdata <- data.frame(wkdata,wk.scalef.obs)
                  glm.Results <- glm(formula(FBoth_two,lhs=2,rhs=2),family=gaussian(link="log"),
                                  data=wkdata,offset=offset.scalef) 
                              } else {
                  glm.Results <- glm(formula(FBoth,lhs=2,rhs=2),family=gaussian(link="log"),
                                  data=wkdata,offset=offset.scalef) } # end if (n.var.0>0) 
               initial.scalef <- coefficients(glm.Results) }
# switching warnings from glm back on
            options(warn=0)                       
            if ((model=='beta binomial') | (model=='correlated binomial')) {
               parameter <- c(initial.p,rep(0,length(glm.Results$coefficients))) 
               } else {
               parameter <- c(initial.p,initial.scalef) 
               } # end of if beta, correlated binomial
            names(parameter) <- c(names(initial.p),names(initial.scalef)) 
            numpar <- length(parameter) } # of if model.type p and scale-factor
                                   } else {
# Checking length of input initial against model value
      parameter <- initial
      numpar    <- length(parameter) 
      if (is.null(names(initial))==TRUE) { 
         cat('\n','WARNING: initial has no associated names','\n')
         names(parameter) <- 1:numpar } # end if is.null(names)
                                      } # end if is.null(initial)

# Checking nobs is the same value as length(list.binary)
   wks <- length(ntrials)
   if (wks!=nobs) { inerror <- 2
         cat('\n','number of rows of covariates.matrix.p not equal to length of list.binary','\n') }
   for ( i in 1:nobs) { nmax1 <- vnmax[i] + 1
                        ntrials[[i]] <- c(rep(vnmax[i],nmax1)) } # end of for loop
   npar.p      <- ncol(covariates.matrix.p)
   npar.scalef <- ncol(covariates.matrix.scalef)
   if (model.type=="p only") {
      npar <- npar.p 
      if (model!="binomial") { npar <- npar + 1 }
      if (numpar!=npar) { inerror <- 2  
            cat('\n','number of parameters error','\n') }}
   if (model.type=="p and scale-factor") {
      npar <- npar.p + npar.scalef
      if (numpar!=npar) { inerror <- 2  
            cat('\n','number of parameters error','\n') }}

# log link function for mean
      r.parameter.p <- rep(0,npar.p) 
      r.parameter.p <- parameter[1:npar.p] 
      lp.p <- covariates.matrix.p%*%r.parameter.p + offset.p
# log-linear function for variance 
      if (model.type=="p and scale-factor") { 
         r.parameter.scalef <- rep(0,npar.scalef) 
         wks <- npar.p + 1
         r.parameter.scalef <- parameter[wks:npar] 
         lp.scalef <- covariates.matrix.scalef%*%r.parameter.scalef + 
                        offset.scalef        } # end of if statement

# start of if inerror if else
if (inerror==0) {   

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
# smallest lower probability allowed 1.e-14
# largest upper probability allowed 1 - 1.e-14
             ll.p <- 1.e-14
             ul.p <- 1 - ll.p
             if (link=='cloglog') { llimit <- log(-log(1-ll.p)) 
                                    ulimit <- log(-log(1-ul.p)) }
             if (link=='logit')   { llimit <- log(ll.p/(1-ll.p)) 
                                    ulimit <- log(ul.p/(1-ul.p)) }
# switching off warnings from optim
             options(warn=-1)
             Results  <- optim(parameter,fn=LL.Regression.Binary,gr=NULL,
                               model.type,model,link,ntrials,nsuccess,
                               covariates.matrix.p,covariates.matrix.scalef,
                               offset.p,offset.scalef,method="Brent",
                               lower=llimit,upper=ulimit,control=control,
                               hessian=FALSE)
# switching warnings from optim back on
             options(warn=0)
                        } else {
              convergence <- list(c("successful"),
                                  c("iteration limit max has been reached"),
                                  c(" "),c(" "),c(" "),c(" "),
                                  c(" "),c(" "),c(" "),c(" "),
                                  c("degeneracy of the Nelder-Mead simplex"))
              Results  <- optim(parameter,fn=LL.Regression.Binary,gr=NULL,
                                model.type,model,link,ntrials,nsuccess,
                                covariates.matrix.p,covariates.matrix.scalef,
                                offset.p,offset.scalef,method="Nelder-Mead",
                                control=control,hessian=FALSE)                            
              wks <- Results$convergence + 1
              cat('\n','optimization method optim:','\n')
              cat(' function calls ',Results$counts[1],'\n')
              cat(' convergence    ',Results$convergence,convergence[[wks]],'\n')
              if (is.null(Results$message)==FALSE) {
                 cat(' message        ',Results$message,'\n') }
                                        } # end of if length(parameter)=1 
          estimates <- Results$par } # end optimization.method=optim
          if (optimization.method=='nlm') {
             options(error=function () { geterrmessage()
                           cat(' Problems calculating the numerical derivative','\n',
                               'nlm is exited. Use optim rather than nlm.','\n') }) # end options 
             Results <- nlm(F.Regression.Binary,parameter,model.type,model,link,ntrials,nsuccess,
                            covariates.matrix.p,covariates.matrix.scalef,
                            offset.p,offset.scalef,
                            hessian=FALSE,fscale=control$fscale,
                            print.level=control$print.level,stepmax=control$stepmax,
                            gradtol=control$gradtol,steptol=control$steptol,
                            iterlim=control$iterlim,ndigit=control$ndigit)
# returning to standard error mode
             options(error=NULL)
             code <- c(
"relative gradient is close to zero,
 current iterate is probably solution",
"succesive iterates within tolerance, 
 current iterate is probably solution",
"last global step failed to locate a point lower than estimate. 
 Either estimate is an approximate lcoal minimum of the function 
 or steptol is too small","iteration limit max has been reached",
"Maximum step size stepmax exceeded five consecutive times. 
 Either the function is unbounded below, becomes asymptotic to a 
 finite value from above in some direction, or stepmax is too small")
             cat('\n','optimization method nlm:','\n')
             cat(' iterations  ',Results$iterations,'\n')                     
             cat(' return code ',Results$code,'\n')
             cat(' ',code[[as.numeric(Results$code)]],'\n') 
             estimates <- Results$estimate } # end optimization.method=nlm 
          names(estimates) <- names(parameter) 
          output <- Model.Binary(estimates,model.type,model,link,ntrials,
                                 covariates.matrix.p,covariates.matrix.scalef,
                                 offset.p,offset.scalef) 
          nobs <- nrow(covariates.matrix.p) 
          mean.par     <- rep(0,nobs) 
          variance.par <- rep(0,nobs) 
          p.par        <- rep(0,nobs)
          scalef.par   <- rep(1,nobs)
          vone         <- rep(1,nobs)
          scalef.limit <- rep(0,nobs)
          exceed.limit <- rep(0,nobs)
# Calculation of p and means from parameter estimates and design matrices
          npar.p      <- ncol(covariates.matrix.p)
          r.parameter.p <- rep(0,npar.p) 
          r.parameter.p <- estimates[1:npar.p] 
# link function for p, logit or cloglog
          lp.p <- covariates.matrix.p%*%r.parameter.p + offset.p
          exp.lp <- exp(lp.p) 
          if (link=='logit')   { p.par <- exp.lp/(vone+exp.lp) }
          if (link=="cloglog") { p.par <- vone - exp(-exp.lp) }

# checking whether the limit for variance of Poisson i.e., variance=mean
# has been reached for any observations for the generalized binomial models
# and that parameters are within limits for generalized, beta and correlated binomials
# obtaining vectors of estimates of theta, rho and limits for beta and correlated binomials
          Dparameters <- output$Dparameters 
# checking that parameters are within limits for generalized, beta and correlated binomials
          if (model.type=="p only") { 
             npar  <- npar.p + 1
             if (model=="generalized binomial") { 
# maximum variance restricted to Poisson i.e., variance=mean
                exceed.limit <- sapply(1:nobs, function(i) 
                             if (Dparameters$vb[i]<=0) { exceed.limit[i] <- 1
                                } else { exceed.limit[i] <- 0 } ) } # end if model
             if ((model=="beta binomial") | (model=="correlated binomial"))
                { theta.llimit <- Dparameters$lower.limit
                  npar <- npar.p + 1
                  vtheta <- rep(estimates[npar],nobs)
                  exceed.limit <- sapply(1:nobs, function(i) 
                            if (vtheta[i]<=theta.llimit[i]) { exceed.limit[i] <- 1
                                                     } else { exceed.limit[i] <- 0 } )
                  if (model=="correlated binomial") { 
                     theta.ulimit <- Dparameters$upper.limit 
                     exceed.limit <- sapply(1:nobs, function(i) 
                            if (vtheta[i]>=theta.ulimit[i]) { exceed.limit[i] <- 1
                                                     } else { exceed.limit[i] <- exceed.limit[i] } ) }
             } # end model=beta or correlated
                                                     } # if model,type

          if (model.type=="p and scale-factor") { 
# Calculation of scale-factors and variances from parameter estimates and design matrices
             npar.scalef <- ncol(covariates.matrix.scalef)
             npar <- npar.p + npar.scalef
# log-linear function for scale-factor 
             r.parameter.scalef <- rep(0,npar.scalef) 
             scalef.par   <- exp(covariates.matrix.scalef%*%r.parameter.scalef + 
                                      offset.scalef)
             if (model=="generalized binomial") { 
# maximum variance restricted to Poisson i.e., variance=mean
                exceed.limit <- sapply(1:nobs, function(i) 
                             if (Dparameters$vb[i]<=0) { exceed.limit[i] <- 1
                                } else { exceed.limit[i] <- 0 } ) } # end if model
             if ((model=="beta binomial") | (model=="correlated binomial"))
                { vtheta <- p.par*(vone-p.par)*scalef.par/(vnmax-vone)
                  theta.llimit <- Dparameters$lower.limit
                  exceed.limit <- sapply(1:nobs, function(i) 
                            if (vtheta[i]<=theta.llimit[i]) { exceed.limit[i] <- 1
                                                     } else { exceed.limit[i] <- 0 } )
                  if (model=="correlated binomial") { 
                     theta.ulimit <- Dparameters$upper.limit 
                     exceed.limit <- sapply(1:nobs, function(i) 
                            if (vtheta[i]>=theta.ulimit[i]) { exceed.limit[i] <- 1
                                                     } else { exceed.limit[i] <- exceed.limit[i] } ) }
                } # end model beta or correlated
                                                 } # if model.type

          ind.exceed.limit <- 0
          if ((model=="generalized binomial") & (sum(exceed.limit)>0)) { ind.exceed.limit <- 1
             if (model.type=="p only") { cat('The parameter b=0 showing that the variance','\n') 
                cat('has reached the Poisson boundary hence its se is set to NA.','\n')
                      } else { cat('The parameter b=0 for some observations showing that','\n') 
                               cat('the variance has reached the Poisson boundary.','\n') }
                             } # end of if model
          if ((model=="beta binomial") & (sum(exceed.limit)>0)) { ind.exceed.limit <- 1
             if (model.type=="p only") { cat('The value of theta is less than the lower limit','\n') 
                cat('for some observations hence its se is set to NA.','\n') 
                      } else { cat('The values of theta for some observations are less','\n') 
                               cat('than the lower limits.','\n') } } # end of if model
          if ((model=="correlated binomial") & (sum(exceed.limit)>0)) { ind.exceed.limit <- 1
             if (model.type=="p only") { cat('The value of rho is outside the lower to upper limit','\n') 
                                         cat('range for some observations hence its se is set to NA.','\n') 
                      } else { cat('The values of rho for some observations are','\n') 
                               cat('outside the lower to upper limit range.','\n') } } # end of if model

# model.hessian from hessian from numDeriv method Richardson
          model.hessian <- hessian(LL.Regression.Binary,
                x=estimates,method="Richardson",
                model.type=model.type,model=model,link=link,
                ntrials=ntrials,nsuccess=nsuccess,
                covariates.matrix.p=covariates.matrix.p,
                covariates.matrix.scalef=covariates.matrix.scalef,
                offset.p=offset.p,offset.scalef=offset.scalef)

            npar <- length(estimates)
# Deleting appropriate row and column of matrix for generalized, beta & correlated binomials
# if limit has been reached
            if (model.type=="p only") { nparm1 <- npar - 1
              if (((model=="generalized binomial") | (model=="beta binomial") |
                   (model=="correlated binomial")) & (ind.exceed.limit>0)) { 
                       nparm1 <- npar - 1
                       nparm1sq <- nparm1*nparm1
                       wk.hessian <- matrix(c(rep(0,nparm1sq)),ncol=nparm1)
                       wk.hessian[nparm1,] <- c(model.hessian[npar,1:nparm1])
                       wk.hessian[,nparm1] <- c(model.hessian[1:nparm1,npar])
                       model.hessian       <- wk.hessian } # if model & ind.exceed.limit
                                                     } # end if model.type
            if (model.type=="p and scale-factor") { 
              if (((model=="generalized binomial") | (model=="beta binomial") |
                   (model=="correlated binomial")) & (ind.exceed.limit>0)) { 
                       npar.psq <- npar.p*npar.p
                       wk.hessian <- matrix(c(rep(0,npar.psq)),ncol=npar.p)
                       model.hessian       <- wk.hessian } # if model & ind.exceed.limit
                                                     } # end if model.type

# checking condition number of model.hessian
          deter <- det(model.hessian)
          if ((is.finite(deter)==FALSE) | (deter==0)) { se <- rep(NA,npar) 
                        } else { listlen <- nrow(covariates.matrix.p)
                                 if ((npar==1) & (listlen==1)) { 
                                    se <- sqrt(-1/model.hessian ) 
                                  } else { condition <- rcond(model.hessian)
# function rcond is from the package Matrix and gives the (reciprocal) condition number
# near 0 is ill-conditioned, near 1 is well conditioned
                                           if (condition>1e-8) { 
                                              varcov <- - solve(model.hessian ) 
                                              se     <- sqrt(diag(varcov)) 
                                                   } else { se <- rep(NA,npar) } } } 

# Replacing appropriate rows of se vector for generalized, beta & correlated binomials
# if limit has been reached
            if (model.type=="p only") {
              if (((model=="generalized binomial") | (model=="beta binomial") |
                   (model=="correlated binomial")) & (ind.exceed.limit>0)) { 
                   se <- c(se[1:nparm1],NA) } # # if model & ind.exceed.limit
                                      } # end if model.type
            if (model.type=="p and scale-factor") {
              if (((model=="generalized binomial") | (model=="beta binomial") |
                   (model=="correlated binomial")) & (ind.exceed.limit>0)) { 
                   se <- c(se[1:npar.p],rep(NA,npar.scalef)) } # # if model & ind.exceed.limit
                                      } # end if model.type

          probabilities <- output$probabilities 
          nobs <- nrow(covariates.matrix.p) 
# Calculation of log likelihood
          loglikelihood <- 0 
          for ( i in 1:nobs) { probability <- probabilities[[i]]
             probability <- probability*(probability>1.e-14) + 1.e-14*(probability<=1.e-14)
             vsuccess <- nsuccess[[i]] 
             if (is.na(probability[1])==TRUE) { loglikelihood <- -1.e+20 
                     } else { loglikelihood  <- loglikelihood + 
                                       t(log(probability))%*%vsuccess }
                             } # end of for loop
     output.fn <- list(model.type=model.type,model=model,link=link,
                       covariates.matrix.p=covariates.matrix.p,
                       covariates.matrix.scalef=covariates.matrix.scalef,
                       offset.p=offset.p,offset.scalef=offset.scalef,
                       estses=data.frame(names(estimates),estimates,se),
                       loglikelihood=loglikelihood,vnmax=vnmax,
                       mean.obs=mean.obs,variance.obs=variance.obs)
                 } else {output.fn <- list(model.type=model.type,model=model,
                           link=link,
                           covariates.matrix.p=covariates.matrix.p,
                           covariates.matrix.scalef=covariates.matrix.scalef,
                           offset.p=offset.p,offset.scalef=offset.scalef,
                           estses=data.frame(NA,NA,NA),loglikelihood=NA,vnmax=NA,
                           mean.obs=mean.obs,variance.obs=variance.obs) } # inerror if else
          return(output.fn)     }
