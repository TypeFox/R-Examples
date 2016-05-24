Model.FaddyJMV.general <-
function(parameter,covariates.matrix.mean,
                    covariates.matrix.variance,offset.mean,offset.variance, 
                    scale.factor.model,fixed.b,vnmax) {
   nobs          <- nrow(covariates.matrix.mean)
   npar.mean     <- ncol(covariates.matrix.mean)
   npar.variance <- ncol(covariates.matrix.variance)
   npar <- npar.mean + npar.variance
   numpar <- length(parameter) 
   r.parameter.mean <- rep(0,npar.mean) 
   r.parameter.mean <- parameter[1:npar.mean] 
   r.parameter.variance <- rep(0,npar.variance) 
   wks <- npar.mean + 1
   r.parameter.variance <- parameter[wks:npar] 
   probabilities <- rep(list(0),nobs)

# log link function for mean
   vlp.mean <- covariates.matrix.mean%*%r.parameter.mean + offset.mean
# log-linear function for variance 
   vlp.variance <- covariates.matrix.variance%*%r.parameter.variance + 
                  offset.variance
   if (scale.factor.model=='yes') { vlp.variance <- vlp.variance + vlp.mean } 
   vlp.mean     <- exp(vlp.mean) 
   vlp.variance <- exp(vlp.variance) 

   if (is.na(fixed.b)==TRUE) {
# extra parameter log(b) 
      wks <- npar + 1
      b   <- exp(parameter[wks])
     logb <- parameter[wks]
                             } else {
# parameter b fixed
      b   <- fixed.b
     logb <- log(b) } # end if.na(fixed.b)

# vectors of Faddy distribution parameters
   out.va <- rep(0,nobs)
   out.vb <- rep(b,nobs)
   out.vc <- rep(0,nobs) 

   for ( i in 1:nobs) {
      nmax  <- vnmax[i]
      nmax1 <- nmax + 1
      vnum  <- c(0:nmax)
      cmean <- vlp.mean[i]
      cvariance <- vlp.variance[i]

# iteration on c
     cmean.b <- cmean+b
     exp.aa  <- cmean.b/b 
     aa <- log(exp.aa) 
     bb <- cvariance/(cmean.b)
     xinc <- ((exp.aa-1-bb)/aa)/((exp.aa*(aa-1)+1)/aa^2)
     if ((((exp.aa-1)>bb)==TRUE) & (is.finite(aa)==TRUE) &
          (is.finite(xinc)==TRUE)) { x <- aa
# limiting case where (exp.aa-1)=bb gives c1=1                              
       old.new <- 0
       iter <- 0
       while ((iter<201) & (abs(xinc)>=1.e-10) & 
              (is.finite(abs(x))==TRUE) & (old.new==0)) {
          old <- x
          x <- x-xinc
          new <- x
          xinc <- ((exp(x)-1)/x-bb/aa)/((exp(x)*(x-1)+1)/x^2)
          iter <- iter + 1
          if (signif(abs(old),14)==signif(abs(new),14)) { old.new=1 }
                             } # end of while
# c <= 1 solution indicated
      c1 <- (x/aa+1)/2 
# rounding c1 to 10 digits  
      c1 <- round(c1,digits=10) 
      if (c1==0) { a1 <- cmean }
      if (c1==1) { a1 <- log(1+cmean/b) }
      if ((c1!=0) & (c1!=1)) { onemc <- 1-c1
          loga1 <- onemc*logb + log(((cmean+b)/b)^onemc-1) - log(onemc)
          if (is.infinite(loga1)==TRUE) { loga1 <- 1.e+308 }
# log(1.e+308) is the last value that returns a number, log(1.e+309) returns Inf                                               
          a1 <- exp(loga1)
#         a1 <- ((cmean+b)^(onemc) - b^(onemc))/onemc
          logva <- rep(loga1,nmax1)
          va <- rep(a1,nmax1)
          vb <- rep(b,nmax1)
          vc <- rep(c1,nmax1)
# vlambda <- va*(vb+vnum)^vc
          vlambda <- loga1 + vc*log(vb+vnum)
          vlambda <- exp(vlambda)
# limiting value for lambda
          lambda.limit <- 745
          vlambda <- sapply(1:nmax1, function(j) 
              if ((is.finite(vlambda[j])==FALSE) | (vlambda[j]>lambda.limit)) { 
                                     vlambda[j] <- lambda.limit  
                                   } else { vlambda[j] <- vlambda[j] } )
                            } # end if ((c1^=0) & (c1^=1))
      out.va[i] <- a1
      out.vc[i] <- c1
# if c1=0 using Poisson function for probabilities
# if c1=1 using negative binomial distribution for probabilities
# only using the EPPM form for c not equal to 0 or 1
      if (c1==0) { probability <- dpois(x=c(0:nmax),lambda=cmean)
                 } else { 
         if (c1==1) { probability <- dnbinom(x=c(0:nmax),size=b,mu=cmean)
                    } else { probability <- EPPMprob(vlambda) } }

          } else { probability <- rep(0,nmax1) } # end of if (((exp.aa-1)...

      probabilities[[i]] <- probability } # end of for loop 
   output <- list(model='general',estimates=parameter,probabilities=probabilities,
                  FDparameters=data.frame(out.va,out.vb,out.vc))
   return(output)              }
