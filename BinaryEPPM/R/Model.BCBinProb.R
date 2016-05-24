Model.BCBinProb <-
function(parameter,model.type,model,link,ntrials,
                        covariates.matrix.p,covariates.matrix.scalef= 
                          matrix(c(rep(1,nrow(covariates.matrix.p))),ncol=1),
                        offset.p=c(rep(0,length(ntrials))),
                        offset.scalef=c(rep(0,length(ntrials)))) {
#  data as number of trials & number of successes  
   twoparameter <- rep(0,2)  
#  the first element of twoparameter is p
#  the second element of twoparameter is theta or rho 
   numpar <- length(parameter) 
   nobs   <- nrow(covariates.matrix.p)
   npar.p      <- ncol(covariates.matrix.p)
   npar.scalef <- ncol(covariates.matrix.scalef)
   npar <- npar.p + npar.scalef
   if (numpar!=numpar) {
      cat('\n','no. of parameters not equal to sum of no of columns mean & variance matrices','\n') }
 
   vone     <- rep(1,nobs) 
   if (model.type=="p and scale-factor") {
      nparm1      <- npar.p + 1
      r.parameter <- rep(0,npar.scalef) 
      r.parameter <- parameter[nparm1:npar]
# modeling the scalefactor with log link
      vscalefact  <- exp(covariates.matrix.scalef%*%r.parameter + offset.scalef) 
                                         } # end if model.type
   r.parameter <- rep(0,npar.p) 
   r.parameter <- parameter[1:npar.p]
   velp <- exp(covariates.matrix.p%*%r.parameter + offset.p) 
   if (link=="logit")   {  vp = velp / (vone + velp) }
   if (link=="cloglog") {  vp = vone - exp(-velp) }
   denom <- rep(0,nobs)
   denom <- sapply(1:nobs, function(i) 
             denom[i] <- max(ntrials[[i]]) )
   vmean     <- denom*vp
   if (model.type=="p and scale-factor") {
      vvariance <- vmean*(vone-vp)*vscalefact
      vtheta <- (vscalefact - vone)/(denom-vone) 
# converting rho to theta for beta binomial so for the
# correlated binomial theta is actually rho
      if (model=="beta binomial") { vtheta <- vtheta/(vone-vtheta) }
      vtheta <- sapply(1:nobs, function(i) 
             if (is.finite(vtheta[i])==FALSE) { vtheta[i] <- 0 
                 } else { vtheta[i] <- vtheta[i] } )
                                         } # end if model.type
   probabilities <- ntrials 
   if (model=="beta binomial") { lower.limit <- rep(0,nobs) }
   if (model=="correlated binomial") { lower.limit <- rep(0,nobs) 
                                       upper.limit <- rep(0,nobs) }
   if (model.type=="p only") { twoparameter[2] <- parameter[npar]
                               vtheta <- rep(twoparameter[2],nobs) } 
   for ( i in 1:nobs) { nt <- denom[i] 
      twoparameter[1] <- vp[i] 
      if (model.type=="p and scale-factor") { twoparameter[2] <- vtheta[i] } 
      if (model=="beta binomial") { 
         wk.output <- BBprob(twoparameter,nt) 
# less than limit for ith observation, setting to limit
         lower.limit[i] <- wk.output$limit } # end of beta binomial
      if (model=="correlated binomial") { 
         wk.output  <- CBprob(twoparameter,nt) 
         lower.limit[i] <- wk.output$limit[1] 
         upper.limit[i] <- wk.output$limit[2] } # end of correlated binomial
      probabilities[[i]] <- wk.output$probability
                      } # end of for loop
      if (model=="beta binomial") { 
          output <- list(model=model,link=link,parameter=parameter,
                         probabilities=probabilities,
                         Dparameters=data.frame(vp,vtheta,lower.limit)) }
      if (model=="correlated binomial") { 
          output <- list(model=model,link=link,parameter=parameter,
                         probabilities=probabilities,
                         Dparameters=data.frame(vp,vrho=vtheta,lower.limit,
                                     upper.limit)) }
   return(output) }
