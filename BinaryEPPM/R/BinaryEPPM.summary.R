BinaryEPPM.summary <-
function(output.fn) {
          options(error=function () { geterrmessage() }) # end options 
# returning to standard error mode
          options(error=NULL)
          cat('\n')
          cat('Model type        :',output.fn$model.type,'\n')
          cat('Model             :',output.fn$model,'\n')
          cat('Link p            :',output.fn$link,'\n')
          if (output.fn$model.type=='p and scale-factor') { 
                         cat('Link scale-factor : log','\n') }
          offsetid.p      <- sum(output.fn$offset.p)
          offsetid.scalef <- sum(output.fn$offset.scalef)
          if ((offsetid.p!=0) | (offsetid.scalef!=0)) { 
              cat('non zero offsets in linear predictors','\n') }   
# checking whether the limit for variance of Poisson i.e., variance=mean
# has been reached for any observations for the generalized binomial models
          nobs <- nrow(output.fn$covariates.matrix.p) 
          mean.par     <- rep(0,nobs) 
          variance.par <- rep(0,nobs) 
          p.par        <- rep(0,nobs)
          scalef.par   <- rep(1,nobs)
          vone         <- rep(1,nobs)
          scalef.limit <- rep(0,nobs)
          exceed.limit <- rep(0,nobs)
# Calculation of p and means from parameter estimates and design matrices
          npar.p      <- ncol(output.fn$covariates.matrix.p)
          r.parameter.p <- rep(0,npar.p) 
          r.parameter.p <- output.fn$estses[,2][1:npar.p] 
# link function for p, logit or cloglog
          lp.p <- output.fn$covariates.matrix.p%*%r.parameter.p + 
                       output.fn$offset.p
          exp.lp <- exp(lp.p) 
          if (output.fn$link=='logit')   { p.par <- exp.lp/(vone+exp.lp) }
          if (output.fn$link=="cloglog") { p.par <- vone - exp(-exp.lp) }

# obtaining vectors of estimates of theta, rho and limits for beta and correlated binomials
          if ((output.fn$model=="beta binomial") | (output.fn$model=="correlated binomial"))
             { ntrials <- list(rep(c(0),nobs))
               for ( i in 1:nobs) { nmax  <- output.fn$vnmax[i]
                                    nmax1 <- nmax + 1 
                             ntrials[[i]] <- c(rep(nmax,nmax1)) } # end of for loop
               output.model <- Model.Binary(parameter=output.fn$estses[,2],
                            model.type=output.fn$model.type,model=output.fn$model,
                            link=output.fn$link,ntrials=ntrials,
                            covariates.matrix.p=output.fn$covariates.matrix.p,
                            covariates.matrix.scalef=output.fn$covariates.matrix.scalef,
                            offset.p=output.fn$offset.p,
                            offset.scalef=output.fn$offset.scalef)
             } # end model

          if (is.na(output.fn$estses[1,3])==TRUE) { 
              cat('Determinant of hessian matrix is zero','\n') 
              cat('or the hessian matrix is ill conditioned.','\n') }

          cat('Parameter estimates and se\'s','\n')
          names(output.fn$estses) <- c('name','Estimates','se')
          print.data.frame(output.fn$estses,row.names=FALSE)
          if (is.na(output.fn$loglikelihood)==TRUE) { cat('loglikelihood is NA','\n')
          } else { cat('\n')
                   cat('log likelihood ',output.fn$loglikelihood,'\n')
                   cat('\n')
                   numpar <- length(output.fn$estses[,2]) 
                   AIC <- -2*output.fn$loglikelihood + 2*numpar 
                   cat('AIC',AIC,'\n')
                 } # end of is.na(output.fn$loglikelihood if
                                   }
