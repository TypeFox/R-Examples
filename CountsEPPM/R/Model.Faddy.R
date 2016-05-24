Model.Faddy <-
function(parameter,model,covariates.matrix.mean,
                offset.mean,fixed.b,vnmax) {
#  data as number of trials & number of successes
#  the parameters c & log(b) are scalars and are the last parameters in the parameter vector   
#  the order of parameters is a, c, b to match with the other functions, b being the nuisance
#  parameter
   nobs <- nrow(covariates.matrix.mean) 
#  covariate matrix for parameter a
   npar.one <- ncol(covariates.matrix.mean)
   npar <- npar.one
   threeparameter <- rep(0,3)  
   va   <- matrix(c(rep(1,nobs)),ncol=1)
   r.parameter <- rep(0,npar.one) 
   r.parameter <- parameter[1:npar.one]
   va <- exp(covariates.matrix.mean%*%r.parameter + offset.mean)
   if (model=="Poisson") { b <- 1
                           c <- 0 } # end of model
   if (model=="negative binomial") { npar <- npar + 1 
                                     b <- exp(parameter[npar])
                                     c <- 1 } # end of model
# parameter b fixed
   if (model=="negative binomial fixed b") { b <- fixed.b 
                                             c <- 1 } # end of model
   if (model=="Faddy distribution") { nparm1 <- npar + 1
                                      npar   <- npar + 2 
                                      c <- parameter[nparm1]
                                      b <- exp(parameter[npar]) } # end of model
# parameter b fixed
   if (model=="Faddy distribution fixed b") { npar <- npar + 1
                                      b <- fixed.b 
                                      c <- parameter[npar] } # end of model
   probabilities <- rep(list(0),nobs)
   threeparameter[2] <- b
   threeparameter[3] <- c 
# vectors of Faddy distribution parameters
   out.va <- va
   out.vb <- rep(b,nobs)
   out.vc <- rep(c,nobs)
   for ( i in 1:nobs) { threeparameter[1] <- va[i] 
                        nmax              <- vnmax[i]
                        nmax1             <- nmax + 1 
                        vid               <- c(0:nmax)
                        if (c<=1) { probability <- Faddyprob.general(threeparameter,nmax)
# c > 1 indicated, all probabilities set to 0
                                  } else { probability <- rep(0,nmax1) } # end of if c<=1
                        probabilities[[i]] <- probability } # end of for i
   output <- list(model=model,estimates=parameter,probabilities=probabilities,
                  FDparameters=data.frame(out.va,out.vb,out.vc))
   return(output)        }
