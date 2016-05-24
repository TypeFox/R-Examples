Model.GB <-
function(parameter,model,link,ntrials,
                        covariates.matrix.p,offset.p=c(rep(0,length(ntrials)))) {
#  data as number of trials & number of successes  
   twoparameter <- rep(0,2)  
#  the second of twoparameter is b 
   npar <- length(parameter) 
   nobs <- nrow(covariates.matrix.p) 
   if (model=="binomial") { twoparameter[2] <- 1 
                            r.parameter     <- rep(0,npar) 
                            r.parameter     <- parameter[1:npar] }
   if (model=="generalized binomial") { nparm1          <- npar - 1 
                                        twoparameter[2] <- parameter[npar] 
                                        r.parameter     <- rep(0,nparm1) 
                                        r.parameter     <- parameter[1:nparm1] }
   vone <- rep(1,nobs) 
   vb   <- rep(twoparameter[2],nobs) 
   velp <- exp(covariates.matrix.p%*%r.parameter + offset.p) 
   if (link=="logit")   { vp = velp / (vone + velp) }
   if (link=="cloglog") { vp = vone - exp(-velp) }
   denom <- rep(0,nobs)
   denom <- sapply(1:nobs, function(i) 
             denom[i] <- max(ntrials[[i]]) )
   vmean <- denom*vp
   probabilities <- ntrials 
   if (round(twoparameter[2],digits=14)==1) { va <- - log(vone - vp)
            } else { vonemb <- vone - vb 
                     va     <- (denom^vonemb - (denom - vmean)^vonemb) / vonemb }
   probabilities <- lapply(1:nobs, function(i) 
          probabilities[[i]] <- GBprob(parameter=c(va[i],vb[i]),nt=denom[i]) )
   output <- list(model=model,link=link,parameter=parameter,
                  probabilities=probabilities,
                  Dparameters=data.frame(va,vb))
   return(output)                                         }
