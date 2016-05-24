BinaryEPPM.distribution <-
function(output.fn,output.probabilities='no',
                                       output.Dparameters='no') {
          options(error=function () { geterrmessage() }) # end options 
# returning to standard error mode
          if (is.na(output.fn$loglikelihood)==FALSE) { 
                   model.type   <- output.fn$model.type
                   model        <- output.fn$model
                   link         <- output.fn$link
                   mean.obs     <- output.fn$mean.obs
                   variance.obs <- output.fn$variance.obs
                   p.obs        <- mean.obs/output.fn$vnmax
                   scalef.obs   <- variance.obs/(mean.obs*(1-p.obs))
                   nobs <- nrow(output.fn$covariates.matrix.p) 
                   ntrials <- list(rep(c(0),nobs))
                   for ( i in 1:nobs) { nmax  <- output.fn$vnmax[i]
                                        nmax1 <- nmax + 1 
                                 ntrials[[i]] <- c(rep(nmax,nmax1)) } # end of for loop
                   mean.par      <- rep(0,nobs) 
                   variance.par  <- rep(0,nobs) 
                   p.par         <- rep(0,nobs)
                   scalef.par    <- rep(1,nobs)
                   mean.prob     <- rep(0,nobs) 
                   variance.prob <- rep(0,nobs) 
                   p.prob        <- rep(0,nobs)
                   scalef.prob   <- rep(1,nobs)
                   vone          <- rep(1,nobs)
                   scalef.limit  <- rep(0,nobs)
                   output.model <- Model.Binary(parameter=output.fn$estses[,2],
                                model.type=model.type,model=model,link=link,
                                ntrials=ntrials,
                                covariates.matrix.p=output.fn$covariates.matrix.p,
                                covariates.matrix.scalef=output.fn$covariates.matrix.scalef,
                                offset.p=output.fn$offset.p,
                                offset.scalef=output.fn$offset.scalef)
                   probabilities <- output.model$probabilities
                   Dparameters   <- output.model$Dparameters

# Calculation of p and means from parameter estimates and design matrices
                   npar.p      <- ncol(output.fn$covariates.matrix.p)
                   r.parameter.p <- rep(0,npar.p) 
                   r.parameter.p <- output.fn$estses[,2][1:npar.p] 
# link function for p, logit or cloglog
                   lp.p <- output.fn$covariates.matrix.p%*%r.parameter.p + 
                              output.fn$offset.p
                   exp.lp <- exp(lp.p) 
                   if (link=='logit')   { p.par <- exp.lp/(vone+exp.lp) }
                   if (link=="cloglog") { p.par <- vone - exp(-exp.lp) }


# checking whether the limit(s), lower for beta, or lower or upper for correlated binomial
# have been reached for any observations 
                   if ((model=="beta binomial") | (model=="correlated binomial")) { 
                      mean.par <- output.fn$vnmax*p.par 
                      scalef.llimit <- vone + Dparameters$lower.limit*(output.fn$vnmax-vone)
                      if (model=="correlated binomial") { 
                         scalef.ulimit <- vone + Dparameters$upper.limit*(output.fn$vnmax-vone) }
                      if (model.type=="p only") { 
                         npar  <- npar.p + 1
                         est.theta  <- output.fn$estses[,2][npar]
                         scalef.par <- vone + rep(est.theta,nobs)*(output.fn$vnmax-vone)
# restricting minimum scale-factor to lower limit
                         scalef.par  <- sapply(1:nobs, function(i) 
                                        if (scalef.par[i]<scalef.llimit[i]) { 
                                            scalef.par[i] <- scalef.llimit[i]
                                        } else { scalef.par[i] <- scalef.par[i] } )
                      if (model=="correlated binomial") { 
                         scalef.par  <- sapply(1:nobs, function(i) 
                                        if (scalef.par[i]>scalef.ulimit[i]) { 
                                            scalef.par[i] <- scalef.ulimit[i]
                                        } else { scalef.par[i] <- scalef.par[i] } ) }
                                                              } # if model,type & model
                      if (model.type=="p and scale-factor") { 
# Calculation of scale-factors and variances from parameter estimates and design matrices
                         npar.scalef <- ncol(output.fn$covariates.matrix.scalef)
                         npar <- npar.p + npar.scalef
# log-linear function for scale-factor 
                         r.parameter.scalef <- rep(0,npar.scalef) 
                         wks <- npar.p + 1
                         r.parameter.scalef <- output.fn$estses[,2][wks:npar] 
                         scalef.par <- exp(output.fn$covariates.matrix.scalef%*%r.parameter.scalef + 
                                                  output.fn$offset.scalef)
# restricting minimum scale-factor to lower limit
                         scalef.par  <- sapply(1:nobs, function(i) 
                                        if (scalef.par[i]<scalef.llimit[i]) { 
                                            scalef.par[i] <- scalef.llimit[i]
                                        } else { scalef.par[i] <- scalef.par[i] } )
# restricting maximum scale-factor to upper limit
                         if (model=="correlated binomial") { 
                            scalef.par  <- sapply(1:nobs, function(i) 
                                        if (scalef.par[i]>scalef.ulimit[i]) { 
                                            scalef.par[i] <- scalef.ulimit[i]
                                        } else { scalef.par[i] <- scalef.par[i] } ) }
                                                             } } # if model=beta binomial

# checking whether the limit for variance of Poisson i.e., variance=mean
# has been reached for any observations for the generalized binomial models
                   if (model=="generalized binomial") { 
                      mean.par <- output.fn$vnmax*p.par 
                      scalef.limit <- vone-p.par
                      if (model.type=="p only") { 
                         npar  <- npar.p + 1
                         wk.2bm1 <- 2*output.fn$estses[,2][npar] - 1
                         scalef.par <- (scalef.limit^wk.2bm1 - 1)/(-p.par*wk.2bm1)
# restricting maximum variance to Poisson i.e., variance=mean
                         scalef.limit <- sapply(1:nobs, function(i) 
                                        if (scalef.limit[i]<1.e-14) { scalef.limit[i] <- 1.e-14
                                          } else { scalef.limit[i] <- scalef.limit[i] } )
                         scalef.par  <- sapply(1:nobs, function(i) 
                                        if (scalef.par[i]>(1/scalef.limit[i])) { 
                                            scalef.par[i] <- 1/scalef.limit[i]
                                        } else { scalef.par[i] <- scalef.par[i] } )
                                                              } # if model,type & model
                      if (model.type=="p and scale-factor") { 
# Calculation of scale-factors and variances from parameter estimates and design matrices
                         npar.scalef <- ncol(output.fn$covariates.matrix.scalef)
                         npar <- npar.p + npar.scalef
# log-linear function for scale-factor 
                         r.parameter.scalef <- rep(0,npar.scalef) 
                         wks <- npar.p + 1
                         r.parameter.scalef <- output.fn$estses[,2][wks:npar] 
                         scalef.par <- exp(output.fn$covariates.matrix.scalef%*%r.parameter.scalef + 
                                                  output.fn$offset.scalef)
# restricting maximum variance to Poisson i.e., variance=mean
                         scalef.limit <- sapply(1:nobs, function(i) 
                                        if (scalef.limit[i]<1.e-14) { scalef.limit[i] <- 1.e-14
                                          } else { scalef.limit[i] <- scalef.limit[i] } )
                         scalef.par  <- sapply(1:nobs, function(i) 
                                        if (scalef.par[i]>(1/scalef.limit[i])) { 
                                            scalef.par[i] <- 1/scalef.limit[i]
                                        } else { scalef.par[i] <- scalef.par[i] } )
                                                             } } # if model=generalized binomial

                   variance.par <- mean.par*(vone-p.par)*scalef.par 

# Calculation of means and variances from count value and predicted probabilities
                      for ( i in 1:nobs) { probability <- probabilities[[i]]
                                    nmax               <- output.fn$vnmax[i] 
                                    vid                <- c(0:nmax)
                                    fmean              <- t(probability)%*%vid 
                                    mean.prob[i]       <- fmean 
                                    variance.prob[i]   <- t(probability)%*%((vid-fmean)^2) 
                                    p.prob[i]          <- mean.prob[i]/nmax
                                    scalef.prob[i]      <- variance.prob[i]/(mean.prob[i]*(1-p.prob[i]))
                                         } # end of for loop

# mean.obs, variance.obs are from the observed means and variances
# mean.par, variance.par are from the parameter estimates and linear predictors
# mean.prob, variance.prob are from the estimated probabilities
                         means.variances <- data.frame(mean.obs,variance.obs,p.obs,scalef.obs,
                                                       mean.par,variance.par,p.par,scalef.par,
                                                       mean.prob,variance.prob,p.prob,scalef.prob)
                   output.distribution <- list(means=means.variances)
                   if (output.probabilities=='yes') { output.distribution <- 
                         c(output.distribution,list(probabilities=probabilities)) } 
                   if (output.Dparameters=='yes') { output.distribution <- 
                         c(output.distribution,list(Dparameters=Dparameters)) } 
                                                      } else { 
                           output.distribution <- list(means=NA,probabilities=NA,
                                                       Dparameters=NA)                                 
                              } # end of is.na(output.fn$loglikelihood if
                   return(output.distribution) }
