Faddyprob.general <-
function(parameter,nmax) {
   nmax1 <- nmax + 1
   vlambda <- rep(0,nmax1) 
   a    <- parameter[[1]] 
   va   <- rep(parameter[[1]],nmax1)
   b    <- parameter[[2]] 
   vb   <- rep(b,nmax1) 
    c   <- parameter[[3]]
   vc   <- rep(c,nmax1)
   vnum <- c(0:nmax)
# c can not be greater than 1 as if it is this function is not
# called in Model.Faddy
# this check is introduced to handle situations where this function
# is called from outside Model.Faddy e.g., when testing the package
# examples
   if (c>1) { cat('\n','c value > 1 as argument to Faddyprob.general',
                  '\n')
              c <- 1 } # end if c>1
   if ((c<1) & (c!=0)) { vlambda <- exp(log(va)+vc*log(vb+vnum))
                                 } # end if c < 1
# if c=0 using Poisson function for probabilities
# if c=1 using negative binomial distribution for probabilities
# only using the EPPM form for c not equal to 0 or 1
      if (c==0) { probability <- dpois(x=vnum,lambda=a)
                 } else { 
         if (c==1) { cmean       <- b*(exp(a) - 1)
                     probability <- dnbinom(x=vnum,size=b,mu=cmean)
                    } else { probability <- EPPMprob(vlambda) } }
   return(probability)   }
