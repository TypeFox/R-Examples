Model.JMVGB <-
function(parameter,model,link,ntrials,
                        covariates.matrix.p,covariates.matrix.scalef,
                        offset.p=c(rep(0,length(ntrials))),
                        offset.scalef=c(rep(0,length(ntrials)))) {
#  data as number of trials & number of successes  
   numpar <- length(parameter) 
   nobs   <- nrow(covariates.matrix.p)
   npar.p      <- ncol(covariates.matrix.p)
   npar.scalef <- ncol(covariates.matrix.scalef)
   npar <- npar.p + npar.scalef
   if (numpar!=numpar) {
      cat('\n','no. of parameters not equal to sum of no. of columns mean & variance matrices','\n') }
   va       <- rep(0,nobs) 
   vb       <- rep(0,nobs) 
   vone     <- rep(1,nobs) 
#  model can only be generalized binomial
   r.parameter <- rep(0,npar.p) 
   r.parameter <- parameter[1:npar.p]
   velp <- exp(covariates.matrix.p%*%r.parameter + offset.p) 
   if (link=="logit")   {  vp = velp / (vone + velp) }
   if (link=="cloglog") {  vp = vone - exp(-velp) }
   denom <- rep(0,nobs)
   denom <- sapply(1:nobs, function(i) 
             denom[i] <- max(ntrials[[i]]) )
   vmean     <- denom*vp
# modeling the scalefactor with log link
   nparm1      <- npar.p + 1
   r.parameter <- rep(0,npar.scalef) 
   r.parameter <- parameter[nparm1:npar]
   vscalefact  <- exp(covariates.matrix.scalef%*%r.parameter + offset.scalef) 
# restricting maximum variance to Poisson i.e., variance=mean
   wkv <- vone-vp
   vscalefact <- sapply(1:nobs, function(i) 
             if (vscalefact[i]>(1/wkv[i])) { vscalefact[i] <- 1/wkv[i]
                } else { vscalefact[i] <- vscalefact[i] } )
   vvariance <- vmean*wkv*vscalefact
   probabilities <- ntrials 
   for ( i in 1:nobs) { nt <- denom[i] 
      vb[i] <- uniroot(function(b,p,scalef) { wk.2bm1 <- 2*b - 1
                        fvalue <- ((1-p)^wk.2bm1 - 1)/(-wk.2bm1*p) - scalef 
                        return(fvalue) },lower=0,upper=1.e+20,
                        p=vp[i],scalef=vscalefact[i],
                        f.lower=(1/(1-vp[i])-vscalefact[i]),
                        f.upper=-vscalefact[i],extendInt="yes",
                        trace=1)$root
      if (round(vb[i],digits=14)==1) { 
                      probabilities[[i]] <- dbinom(c(0:denom[i]),denom[i],vp[i],log=FALSE)
             } else { onemb <- 1 - vb[i]
                      va[i] <- (denom[i]^onemb - (denom[i] - vmean[i])^onemb) / onemb
                      wk.prob <- GBprob(parameter=c(va[i],vb[i]),nt=denom[i]) 
                      probabilities[[i]] <- wk.prob }
                    } # end of for loop
   output <- list(model=model,link=link,parameter=parameter,
                  probabilities=probabilities,
                  Dparameters=data.frame(va,vb))
   return(output)                                         }
