# This function returns the solution of the minimization problem
setGeneric("zsuml2min",
           function (o, x=0, y=0, max.iter=100, eps=1.e-3, verbose=FALSE, algorithm="weiszfeld", ...) standardGeneric("zsuml2min")
)

# General zsuml2min function
# L-BFGS-B seems to be the best similar to weiszfeld
# Take into account that weiszfeld is completely implemented in R
setMethod("zsuml2min", "loca.p",
function (o, x=0, y=0, max.iter=100, eps=1.e-3, verbose=FALSE, algorithm="weiszfeld", ...)
   {
   if (algorithm=="gradient" || algorithm=="g") zsuml2mingradient.loca.p(o, x, y, max.iter, eps, verbose)
   else if (algorithm=="search" || algorithm=="s") zsuml2minsearch.loca.p(o, x, y, max.iter, eps, verbose)
   else if (algorithm=="weiszfeld" || algorithm=="w") zsuml2minweiszfeld.loca.p(o, x, y, max.iter, eps, verbose, ...)
   else if (algorithm=="ucminf" || algorithm=="u") zsuml2minucminf.loca.p(o, x, y, max.iter, eps, verbose)
   else
     {
       zzsummin <- function(x) zsum(o, x[1], x[2])
       par <- c(sum(o@x*o@w)/sum(o@w), sum(o@y*o@w)/sum(o@w)) 
       optim(par, zzsummin, method=algorithm, control=list(maxit=max.iter))$par
     }
   }
)

# Optimization by ucminf function from ucminf package
zsuml2minucminf.loca.p <- function (o, x=0, y=0, max.iter=100, eps=1.e-3, verbose=FALSE)
   {
     require('ucminf')
     zzsum <- function(xx) zsum(o, xx[1], xx[2])
     sol <- ucminf(par = c(x, y), fn = zzsum, control=list(maxeval=max.iter, trace=verbose))
     if (verbose) cat(gettext(sol$message));
     return(sol$par)
   }

# Gradient Method
zsuml2mingradient.loca.p <- function (o, x=0, y=0, max.iter=100, eps=1.e-3, verbose=FALSE)
   {
   lambda <- 1;
   eps2 <- eps^2
   u<-c(x,y)
   z <- zsum(o, u[1], u[2])
   for (i in 0:max.iter)
      {
      if (verbose) cat(paste(gettext("Iter.", domain = "R-orloca"), i, ": (", u[1], ",", u[2], ") ", z, "\n", sep=""))
      g<-zsumgra(o, u[1], u[2])
      mg <- sum(g^2)
      # Check stop rule
      if (is.na(mg))
        {
        # A demand point stop rule
        g<-zsumgra(o, u[1], u[2], partial=T)
        mg <- sum(g^2)
        ii <- which.min((o@x-u[1])^2+(o@y-u[2])^2)
        if (mg < sum(o@w[ii]^2)) 
        	  {
        	  if(verbose) cat(gettext("Optimality condition reached at demand point.", domain = "R-orloca"));
        	  break
        	  }
        }
      else if (mg < eps2)
         {
         if(verbose) cat(gettext("Optimality condition reached.", domain = "R-orloca"));
         break;
         }
      nu <- u - lambda*g
      nz <- zsum(o, nu[1], nu[2])
      if (nz < z)
         {
         u<-nu
         z<-nz
         lambda <- lambda*2.2
         }
      else
         {
         lambda <- lambda/2
         }
      }
   if (i == max.iter) warning.max.iter(max.iter)
   u
   }

zsuml2minsearch.loca.p <- function (o, x=0, y=0, max.iter=100, eps=1.e-3, verbose=FALSE)
   {
   warning(gettext('Deprecated option for algorithm in zsummin.', domain = "R-orloca"))
   eps2 <- eps^2
   lambda <- c(1, 1)
   u <- c(x, y)
   z <- zsum(o, x, y)
   nu <- u
   for(i in 0:max.iter)
      {
      for (j in 1:2)
         {
         nu[j] <- u[j] + lambda[j]
         nz <- zsum(o, nu[1], nu[2])
         if (nz < z)
            {
            u <- nu
            z <- nz
            lambda[j] <- 2.2 * lambda[j]
            }
         else
            {
            nu[j] <- u[j] - lambda[j]
            nz <- zsum(o, nu[1], nu[2])
            if (nz < z)
               {
               u <- nu
               z <- nz
               lambda[j] <- -2.2 * lambda[j]
               }
            else lambda[j] <- lambda[j]/2
            }

         }
      if (verbose) cat(paste(gettext("Iter.", domain = "R-orloca"), i, ": (", u[1], ",", u[2], ") ", z, "\n", sep=""))
      if (sum(lambda^2) < eps2)
        {
        if(verbose) cat(gettext("Optimality condition reached.", domain = "R-orloca"));
        break;
        }
      }
   if (i == max.iter) warning.max.iter(max.iter)
   u
   }

zsuml2minweiszfeld.loca.p <- function (o, x=0, y=0, max.iter=100, eps=1.e-3, verbose=FALSE, csmooth=.9)
   {
   # Check smooth value
   if (!identical(csmooth >= 0 && csmooth < 1, TRUE))
     {
       warning(paste(gettext("Value for smooth parameter non valid:", domain = "R-orloca"), smooth, gettext("Reseting to its default value.", domain = "R-orloca")))
       csmooth <- .5
     }
   eps2 <- eps^2
   u<-c(x,y)
   # Begin iterations in non smooth mode
   .smooth = 0
   i.i = 0
   i.s = round(max.iter*.5)
   for (j in 1:2)
     {
   for (i in i.i:i.s)
      {
      if (verbose) cat(paste(gettext("Iter. ", domain = "R-orloca"), i, ": (", u[1], ",", u[2], ") ", zsum(o, u[1], u[2]), "\n", sep=""))
      # Compute the distances to demand points
      n <- sqrt((u[1]-o@x)^2+(u[2]-o@y)^2)
      # Check for demand point proximities
      ii <- (n > eps)
      # Compute the numerator of iteration
      n <- o@w/n;
      # Compute the gradient
      g <- c(sum((u[1]-o@x[ii])*n[ii]), sum((u[2]-o@y[ii])*n[ii]))
      mg <- sum(g^2)
     # Check stop rule
      if (all(ii))
         {
         # A demand point stop rule
         if (mg < sum(o@w[!ii]^2) || mg < eps2)
           {
           if(verbose) cat(gettext("Optimality condition reached at demand point.", domain = "R-orloca"));
           break
           }
         }
      # Generic stop rule
      else if (mg <eps2)
        {
        if(verbose) cat(gettext("Optimality condition reached.", domain = "R-orloca"));
        break
        }
      
      s <- sum(n[ii])
      nx <- n*o@x
      ny <- n*o@y
      u <- .smooth * u + (1-.smooth) * c(sum(nx[ii]), sum(ny[ii]))/s
      }
      # Check if optimality condition had been reached
      if (i != i.s) break
      # Changing to smooth version
      .smooth = csmooth
      if (j == 1) warning(gettext("The algorithm seems converges very slowly. Trying now with the smooth version.", domain = "R-orloca"))
   i.i = i.s
   i.s = max.iter
 }
   if (i == max.iter) warning.max.iter(max.iter)
   u
   }
