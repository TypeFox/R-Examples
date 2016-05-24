# This function returns the solution of the minimization problem
setGeneric("zsumlpmin",
           function (o, x=0, y=0, p=2, max.iter=100, eps=1.e-3, verbose=FALSE, algorithm="weiszfeld", ...) standardGeneric("zsumlpmin")
)

# Optimization by ucminf function from ucminf package
zsumlpminucminf.loca.p <- function (o, x=0, y=0, p=2, max.iter=100, eps=1.e-3, verbose=FALSE)
   {
     require('ucminf')
     zzsum <- function(xx) zsumlp(o, xx[1], xx[2], p=p)
     sol <- ucminf(par = c(x, y), fn = zzsum, control=list(maxeval=max.iter, trace=verbose))
     if (verbose) cat(gettext(sol$message));
     return(sol$par)
   }

# Gradient Method
setMethod("zsumlpmin", "loca.p",
function (o, x=0, y=0, p=2, max.iter=100, eps=1.e-3, verbose=FALSE, algorithm="weiszfeld", ...)
   {
     if (p>=1) {
       if (algorithm=="gradient" || algorithm=="g") zsumlpmingradient.loca.p(o, x, y, p, max.iter, eps, verbose)
       else if (algorithm=="search" || algorithm=="s") zsumlpminsearch.loca.p(o, x, y, p, max.iter, eps, verbose)
       else if (algorithm=="weiszfeld" || algorithm=="w") zsumlpminweiszfeld.loca.p(o, x, y, p, max.iter, eps, verbose, ...)
       else if (algorithm=="ucminf" || algorithm=="u") zsumlpminucminf.loca.p(o, x, y, p, max.iter, eps, verbose)
   else
     {
       zzsummin <- function(x) zsumlp(o, x[1], x[2], p=p)
       par <- c(sum(o@x*o@w)/sum(o@w), sum(o@y*o@w)/sum(o@w)) 
       optim(par, zzsummin, method=algorithm, control=list(maxit=max.iter))$par
     }
     }
   else stop(paste(p, gettext("is not a valid value for p, use 1 <= p", domain = "R-orloca")))
   }
)

zsumlpmingradient.loca.p <- function (o, x=0, y=0, p=2, max.iter=100, eps=1.e-3, verbose=FALSE)
   {
   lambda = 1;
   eps2 <- eps^2
   u<-c(x,y)
   z <- zsumlp(o, u[1], u[2], p)
   for (i in 0:max.iter)
      {
      if (verbose) cat(paste(gettext("Iter.", domain = "R-orloca"), i, ": (", u[1], ",", u[2], ") ", z, "\n", sep=""))
      g <- zsumlpgra(o, u[1], u[2], p)
      mg <- sum(g^2)
      if (is.na(mg))
         {
         # A demand point stop rule
         g <- zsumlpgra(o, u[1], u[2], p, partial=T)
         q <- p/(p-1)
         mg <- sum(abs(g)^q)^(1/q)
         ii <- which.min((o@x-u[1])^2+(o@y-u[2])^2)
         if (mg < sum(o@w[ii]))
        	  {
        	  if(verbose) cat(gettext("Optimality condition reached at demand point.", domain = "R-orloca"));
        	  break
        	  }
         }
      else if (mg<eps2)
         {
         if(verbose) cat(gettext("Optimality condition reached.", domain = "R-orloca"));   
         break;
         }
      nu <- u - lambda*g
      nz <- zsumlp(o, nu[1], nu[2], p)
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
   if (verbose && i == max.iter) cat(gettext("Maximun number of iteration reached"));
   u
   }

zsumlpminsearch.loca.p <- function (o, x=0, y=0, p=2, max.iter=100, eps=1.e-3, verbose=FALSE)
   {
   warning(gettext('Deprecated option for algorithm in zsummin', domain = "R-orloca"))
   eps2 <- eps^2
   lambda <- c(1, 1)
   u <- c(x, y)
   z <- zsumlp(o, x, y, p)
   nu <- u
   for(i in 0:max.iter)
      {
      for (j in 1:2)
         {
         nu[j] <- u[j] + lambda[j]
         nz <- zsumlp(o, nu[1], nu[2], p)
         if (nz < z)
            {
            u <- nu
            z <- nz
            lambda[j] <- 2.2 * lambda[j]
            }
         else
            {
            nu[j] <- u[j] - lambda[j]
            nz <- zsumlp(o, nu[1], nu[2], p)
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
        break
        }
      }
   if (verbose && i == max.iter) cat(gettext("Maximun number of iteration reached", domain = "R-orloca"));
   u
   }

zsumlpminweiszfeld.loca.p <- function (o, x=0, y=0, p=2, max.iter=100, eps=1.e-3, verbose=FALSE, csmooth=.5)
   {
   # Check smooth value
   if (!identical(csmooth >= 0 && csmooth < 1, TRUE))
     {
       warning(paste(gettext("Value for smooth parameter non valid", domain = "R-orloca"), smooth, gettext("Reseting to its default value", domain = "R-orloca")))
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
      if (verbose) cat(paste(gettext("Iter. ", domain = "R-orloca"), i, ": (", u[1], ",", u[2], ") ", zsumlp(o, u[1], u[2], p), "\n", sep=""))
      # Compute the distances to demand points
      n <- (abs(u[1]-o@x)^p+abs(u[2]-o@y)^p)^(1/p)
      # Check for demand point proximities
      ii <- (n > eps)
      # Compute the numerator of iteration
      n <- o@w*(abs(u[1]-o@x)^p+abs(u[2]-o@y)^p)^(1/p-1)
      # Compute the gradient
      g <- c(sum(sign(u[1]-o@x[ii])*abs(u[1]-o@x[ii])^(p-1)*n[ii]), sum(sign(u[2]-o@y[ii])*abs(u[2]-o@y[ii])^(p-1)*n[ii]))
      mg <- sum(g^2)
      # Check stop rule
      if (all(ii))
         {
         # A demand point stop rule
         q <- p/(p-1)
         mg <- sum(abs(g)^q)^(1/q)
         if (mg < sum(o@w[!ii]) || mg < eps2)
           {
           if(verbose) cat(gettext("Optimality condition reached at demand point.", domain = "R-orloca"));
           break
           }
         }
      # Generic stop rule
      else if (mg<eps2)
        {
        if(verbose) cat(gettext("Optimality condition reached.", domain = "R-orloca"));
        break;
        }
      dx <- n*abs(u[1]-o@x)^(p-2)
      nx <- dx*o@x
      dy <- n*abs(u[2]-o@y)^(p-2)
      ny <- dy*o@y
      u <- .smooth * u + (1-.smooth) * c(sum(nx[ii])/sum(dx[ii]), sum(ny[ii])/sum(dy[ii]))
      }
      # Check if optimality condition had been reached
      if (i != i.s) break
      # Changing to smooth version
      .smooth = csmooth
      if (j == 1) warning(gettext("The algorithm seems converges very slowly. Trying now with the smooth version."))
   i.i = i.s
   i.s = max.iter

 }
   if (i == max.iter) warning.max.iter(max.iter)
   u
   }
