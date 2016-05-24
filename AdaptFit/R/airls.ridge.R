###### R function: airls.ridge ##########

# For fitting an iteratively reweighted
# least squares ridge regression.
# The C=[X Z] notation for mixed model
# representation is used in this function.


airls.ridge <- function(C.mat,y,off.var=rep(0,length(y)),ridge.vec,
                       max.it,acc,family="binomial",track=FALSE)
{
   # First fit an ordinary ridge regression.
   # to obtain initial values.

   y.off <- y - off.var

   ridge.reg.fit <- ridge.reg(C.mat,y.off,ridge.vec=ridge.vec)
   beta.hat <- ridge.reg.fit$coef

   muhat <- off.var + C.mat%*%beta.hat

   muhat <-  gen.corr.range(muhat,family)
 
   etahat <- link(muhat,family)

   # Now apply iteratively reweighted least squares
   # to fit generalized version.

   yhat.old <- inv.link(etahat,family)

   converged <- FALSE
   it.num <- 0
   while (converged==FALSE) 
   {
      it.num <- it.num + 1
      if (it.num>max.it) stop("maximum number of iterations exceeded")
   
      # Obtain the weights
   
      muhat <- inv.link(etahat,family)
      wt <- dinv.link(etahat,family)

      # Obtain the adjusted dependent variable
   
      z <- etahat + (y - muhat)/wt  
   
      # Update the fit
   
      z.off <- z - off.var 
   
      ridge.reg.fit <- ridge.reg(C.mat,z.off,wt,ridge.vec=ridge.vec)
      beta.hat <-  ridge.reg.fit$coef
   
      etahat <- off.var + C.mat%*%beta.hat 
   
      # Check for convergence
   
      yhat.new <- inv.link(etahat,family)
 
      err <- sqrt(sum((yhat.new-yhat.old)^2))/sqrt(sum(yhat.old^2))

      if (err<acc) converged <- TRUE
   
      if (track==TRUE)
      {
         cat("\n")
         cat("Current error in IRLS is",round(err,8),"\n")
         cat("\n")  
      }
   
      yhat.old <- yhat.new
   }
   
   return(ridge.reg.fit)
}

