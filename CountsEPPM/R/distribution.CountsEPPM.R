distribution.CountsEPPM <-
function(output.fn,output.probabilities='no',
                                       output.FDparameters='no') {
          if (is.finite(output.fn$loglikelihood)==TRUE) { 
                    output <- Model.Counts(parameter=output.fn$estses[,2],
                              model.type=output.fn$model.type,
                              model=output.fn$model,
                              covariates.matrix.mean=output.fn$covariates.matrix.mean,
                              covariates.matrix.variance=output.fn$covariates.matrix.variance,
                              offset.mean=output.fn$offset.mean,
                              offset.variance=output.fn$offset.variance,
                              scale.factor.model=output.fn$scale.factor.model,
                              fixed.b=output.fn$fixed.b,vnmax=output.fn$vnmax) 
                   probabilities <- output$probabilities
                   FDparameters  <- output$FDparameters
                   nobs <- nrow(output.fn$covariates.matrix.mean) 
# Truncation
                   ltvalue <- output.fn$ltvalue
                   utvalue <- output.fn$utvalue
                   for ( i in 1:nobs) { probability <- probabilities[[i]]     
                      nmax1 <- output.fn$vnmax[i] + 1
                      if ((is.na(ltvalue)==FALSE) | 
                          (is.na(utvalue)==FALSE)) { 
                          rev.probability <- LRTruncation(probability,ltvalue,utvalue)
                          if (is.na(ltvalue)==FALSE) { wks1 <- ltvalue + 2 
                                                     } else { wks1 <- 1 }
                          if (is.na(utvalue)==FALSE) { wks2 <- utvalue  
                                                     } else { wks2 <- nmax1 }
                          probabilities[[i]] <- rep(0,nmax1)
                          probabilities[[i]][wks1:wks2] <- rev.probability[wks1:wks2] }
                                       } # end of for loop
                      mean.obs       <- output.fn$mean.obs
                      variance.obs   <- output.fn$variance.obs
                      scalef.obs     <- rep(0,nobs) 
                      mean.par       <- rep(0,nobs) 
                      variance.par   <- rep(0,nobs) 
                      scalef.par     <- rep(0,nobs) 
                      mean.prob      <- rep(0,nobs) 
                      variance.prob  <- rep(0,nobs) 
                      scalef.prob    <- rep(0,nobs) 
                      totalprob      <- rep(0,nobs)
                      if (output.fn$model.type=="mean and variance") { 
# Calculation of means and variances from parameter estimates and design matrices
                      npar.mean     <- ncol(output.fn$covariates.matrix.mean)
                      npar.variance <- ncol(output.fn$covariates.matrix.variance)
                      npar <- npar.mean + npar.variance
                      r.parameter.mean <- rep(0,npar.mean) 
                      r.parameter.mean <- output.fn$estses[,2][1:npar.mean] 
# log link function for mean
                      lp.mean <- output.fn$covariates.matrix.mean%*%r.parameter.mean + 
                                 output.fn$offset.mean
                      mean.par <- exp(lp.mean)
# log-linear function for variance 
                      r.parameter.variance <- rep(0,npar.variance) 
                      wks <- npar.mean + 1
                      r.parameter.variance <- output.fn$estses[,2][wks:npar] 
                      lp.variance <- output.fn$covariates.matrix.variance%*%r.parameter.variance + 
                                     output.fn$offset.variance 
                      if (output.fn$scale.factor.model=='yes') { lp.variance <- 
                                   lp.variance + log(mean.par) }
                      variance.par <- exp(lp.variance) } # if model.type
# Calculation of means and variances from count value and predicted probabilities
                      for ( i in 1:nobs) { probability <- probabilities[[i]]
                                    nmax              <- output.fn$vnmax[i] 
                                    vid               <- c(0:nmax)
                                    fmean             <- t(probability)%*%vid 
                                    mean.prob[i]      <- fmean 
                                    variance.prob[i]  <- t(probability)%*%((vid-fmean)^2) 
                                    totalprob[i] <- sum(probability) } # end of for loop
                      if (output.fn$model.type=="mean and variance") {
# mean.obs, variance.obs, scalef.obs are the observed means & variances from the data
# mean.par, variance.par, scalef.par are from the parameter estimates & linear predictors
# mean.prob, variance.prob, scalef.prob are from the estimated probabilities
                         scalef.obs  <- variance.obs/mean.obs
                         scalef.par  <- variance.par/mean.par
                         scalef.prob <- variance.prob/mean.prob
                         means.variances <- data.frame(mean.obs,variance.obs,scalef.obs,
                                              mean.par,variance.par,scalef.par,
                                              mean.prob,variance.prob,scalef.prob,totalprob) 
                         } else { 
                         means.variances <- data.frame(mean.obs,variance.obs,
                                                       mean.prob,variance.prob,totalprob)
                                 } # end if model.type
                      if (output.probabilities=='yes') { 
                         output.distribution <- list(means=means.variances,probabilities=probabilities)
                          } else { 
                         output.distribution <- list(means=means.variances) } # end if
                      if (output.FDparameters=='yes') { 
                         output.distribution <- c(output.distribution,
                                                  list(FDparameters=FDparameters)) } # end if
                                                         } else { 
                          output.distribution <- list(means=NA,probabilities=NA,FDparameters=NA)                                 
                                 } # end of is.na(output.fn$loglikelihood if
                      return(output.distribution) }
