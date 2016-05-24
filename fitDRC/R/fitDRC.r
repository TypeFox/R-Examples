################################################################################
# Required elements of class "distribution":
# ------------------------------------------
#
# name            character        name
# range           function         range over which the distribution is defined
#                                  (bounds can be -Inf or +Inf)
# par.names       character(npar)  parameter names (npar names)
# par.ranges      matrix(npar,2)   ranges of parameters (one after the other)
# par             numeric(npar)    (default) parameter values (npar values)
# pdf             function         pdf(x,par)
# cdf             function         cdf(x,par)
# cdf.inv         function         cdf.inv(p,par)
# mean            function         mean(par)
# sd              function         sd(par)
# mode            function         mode(par)
# median          function         median(par)
#
# Methods of class "distribution":
# --------------------------------
#
# print(dist,...)
# summary(dist,...)
# plot(dist,par, range=NA, plot=TRUE,...)     produces a pdf listing if plot=FALSE

print.distribution <- function(x = dist,...)
{  dist <- x
   rm(x)
   title <- paste(dist$name,"Distribution")
   cat(paste(title,"\n"))
   cat(paste(paste(rep("-",nchar(title)),collapse=""),"\n"))
   cat(paste("Range:",RANGE(dist)[1],RANGE(dist)[2],"\n"))
   cat("Parameter names, ranges and values:\n")
   ranges <- dist$par.ranges
   tab <- cbind(ranges,dist$par)
   colnames(tab) <- c("Min","Max","Value")
   rownames(tab) <- dist$par.names
   print(tab)
   cat("\n")
   cat(paste("Mean:  ",MEAN(  dist),"\n"))
   cat(paste("StDev: ",SD(    dist),"\n"))
   cat(paste("Median:",MEDIAN(dist),"\n"))
   cat(paste("Mode:  ",MODE(  dist),"\n"))
   
}

summary.distribution <- function(object = dist,...)
{
   print(object,...)
}

plot.distribution <- function(x = dist, par=dist$par, range=NA, what="PDF", plot=TRUE, length=101,...)
{  dist <- x
   rm(x)
   if ( what == "PDF" | what == "CDF" )
   {
      if ( is.na(range[1]) )
      {
         x.min <- CDFinv(dist,0.001,par)
         x.max <- CDFinv(dist,0.999,par)
      }
      else
      {
         x.min <- range[1]
         x.max <- range[2]
      }
      x <- seq( x.min, x.max, length = length )
      if ( what == "PDF" )
      {
         y <- PDF(dist,x,par)
      }
      else
      {
         y <- CDF(dist,x,par)
      }
   }
   else
   {
      if ( what == "CDFinv" )
      {
         if ( is.na( range[1]) )
         {
            x.min <- 0.001
            x.max <- 0.999
         }
         else
         {
            x.min <- range[1]
            x.max <- range[2]
         }
         x <- seq( x.min, x.max, length = length )
         y <- CDFinv( dist, x, par )
      }
      else
      {
         stop("*** illegal value of what in plot.distribution")
      }
   }   
   if ( plot )
   {
      plot( x, y, type="l",... )
   }
   else
   {
      res <- cbind(x,y)
      colnames(res) <- c("x","y")
      return(res)
   }
}

mergePar <- function(par, par.default, par.names=NA,...)
{
   # check consistency of input; 
   # initialize parameter vector p with default values:
   
   n.default <- par.names
   if ( is.na(par.names[1]) ) n.default <- names(par.default)
   if ( length(par.default) != length(n.default) ) 
   {
      stop("*** default parameter values and names not of the same length")
   }
   
   p <- par.default
   names(p) <- n.default

   # replace selected default values with values provided by par
   
   n <- names(par)
   if ( length(n) == 0 ) # no labels required if of full length (correct order assumed)
   {
      if ( length(par) == length(p) )
      {
         p <- par
         names(p) <- n.default
      }
      else
      {
         stop("*** provided parameter vector must be labelled if not of full length")
      }
   }
   else
   {
      for ( i in 1:length(n) ) # replace values with given names by provided values
      {
         m <- match(n[i],n.default)
         if ( !is.na(m) )
         {
            p[m] <- par[i]
         }
      }
   }
   return(p)
}

PDF.distribution <- function(dist, x, par = dist$par, log = FALSE,...)
{
   par.local <- mergePar(par, dist$par, dist$par.names)
   
   # if no transformation, just calculate pdf:
   if ( length(dist$trans) == 0 )
   {
      if ( !log )
      {
      	return(dist$pdf(x,par.local))
      }
      else # return the logarithm
      {
      	return( log(dist$pdf(x,par.local))  )
      }
   	}
   # with transformation calculate f_x(x) = f_t(t(x)) * dt/dx(x) :
   else
   {
   npar.dist <- length( dist$dist$par ) 
   par.dist  <- par.local[  1:npar.dist ]
   par.trans <- par.local[-(1:npar.dist)] 
      if ( !log )
      { 
   return(dist$dist$pdf(dist$trans$trans.forw(x,par.trans),par.dist) *
          dist$trans$trans.deriv( x,par.trans ) )
       }
       else # return the logarithm
       {
       return( log(dist$dist$pdf(dist$trans$trans.forw(x,par.trans),par.dist) *
               dist$trans$trans.deriv( x,par.trans )) )
       }
   }
}

CDF.distribution <- function(dist, x, par = dist$par, ...)
{
   par.local <- mergePar( par, dist$par, dist$par.names )
   
   # if no transformation, just calculate cdf:
   if ( length(dist$trans) == 0 ) return(dist$cdf(x,par.local))   

   # with transformation calculate F_x(x) = F_t(t(x)) :

   npar.dist <- length(dist$dist$par) 
   par.dist  <- par.local[  1:npar.dist ]
   par.trans <- par.local[-(1:npar.dist)]
   return(dist$dist$cdf(dist$trans$trans.forw(x,par.trans),par.dist))
}

CDFinv.distribution <- function(dist, p, par = dist$par,...)
{
   par.local <- mergePar(par,dist$par,dist$par.names)

   # if no transformation, just calculate cdf.inv:
   if ( length(dist$trans) == 0 ) return(dist$cdf.inv(p,par.local))
   
   # with transformation calculate F_x.inv = t.inv(F_t.inv(p)) :
   npar.dist <- length(dist$dist$par) 
   par.dist  <- par.local[  1:npar.dist ]
   par.trans <- par.local[-(1:npar.dist)]
   return(dist$trans$trans.backw(dist$dist$cdf.inv(p,par.dist),par.trans))
}


RANGE.distribution <- function(dist, par = dist$par,...)
{
   par.local <- mergePar(par,dist$par,dist$par.names)

   # if no transformation, just return the range of the distribution:
   if ( length(dist$trans) == 0 )
   {
      return(dist$range(par.local))
   }
   # with transformation return the ranges of the transformation:
   else
   {
      npar.dist <- length(dist$dist$par) 
      par.dist  <- par.local[  1:npar.dist ]
      par.trans <- par.local[-(1:npar.dist)]
      return(dist$trans$range.x(par.trans))
   }
}


MEAN.distribution <- function(dist, par = dist$par,...)
{
	par.local <- mergePar(par,dist$par,dist$par.names)
	ifelse(
		is.null(dist$mean),
		return(NA),
		return(dist$mean(par.local))
	)
}

SD.distribution <- function(dist, par = dist$par,...)
{
	par.local <- mergePar(par, dist$par, dist$par.names)
	ifelse(
		is.null(dist$sd),
		return(NA),
		return(dist$sd(par.local))
	)
}

MEDIAN.distribution <- function(dist, par = dist$par,...)
{
	par.local <- mergePar(par,dist$par,dist$par.names)
	ifelse(
		is.null(dist$median),
		return(NA),
		return(dist$median(par.local))
	)
}

MODE.distribution <- function(dist, par = dist$par,...)
{
	par.local <- mergePar(par,dist$par,dist$par.names)
	ifelse(
		is.null(dist$mode),
		return(NA),
		return(dist$mode(par.local))
	)
}

PDF     <- function(dist,x,par,...) UseMethod("PDF")
CDF     <- function(dist,x,par,...) UseMethod("CDF")
CDFinv  <- function(dist,p,par,...) UseMethod("CDFinv")
RANGE   <- function(dist,par,  ...) UseMethod("RANGE")
MEAN    <- function(dist,par,  ...) UseMethod("MEAN")
SD      <- function(dist,par,  ...) UseMethod("SD")
MEDIAN  <- function(dist,par,  ...) UseMethod("MEDIAN")
MODE    <- function(dist,par,  ...) UseMethod("MODE")

################################################################################

# Normal distribution constructor:
# --------------------------------

dist.normal.create <- function(par=NA)
{
   # set default parameter values:

   par.default <- c(0,1)
   names(par.default)     <- c("Mean","StDev")

   # overwrite parameter values with provided values:

   p <- mergePar(par,par.default)

   # construct class:
   
   dist            <- list()
   dist$name       <- "Normal"
   dist$range      <- function(par)
                      {
                         return(c(-Inf,+Inf))
                      }
   dist$par.names  <- names(p)
   dist$par.ranges <- matrix(c(-Inf,+Inf,  # Mean
                               0,Inf),     # Standard Deviation
                             byrow=TRUE,ncol=2)
   dist$par        <- p
   dist$mean       <- function(par)
                      {
   	                   mean <- par[1]
   	                   names(mean) <- "Mean"
   	                   return(mean)
   	                   }
   dist$sd         <- function(par)
                      {
                         sd <- par[2]
                         names(sd) <- "StDev"
                         return(sd)
                      }
   dist$median     <- function(par)
                      {
   	                   median <- par[1]
   	                   names(median) <- "Median"
   	                   return(median)
   	                   }
   dist$mode       <- function(par)
                      {
                      mode <- par[1]
                      names(mode) <- "Mode"
                      return(mode)
                      }
   dist$pdf        <- function(x,par) { return(dnorm(x,par[1],par[2])) }
   dist$cdf        <- function(x,par) { return(pnorm(x,par[1],par[2])) }
   dist$cdf.inv    <- function(p,par) { return(qnorm(p,par[1],par[2])) }
   class(dist)     <- "distribution"

   return(dist)
}

################################################################################

# Tests:
# ------
# print(dist.normal.create(c(Mean = 0, StDev = 1)))
# print(dist.normal.create(c(Mean  = 99)))
# print(dist.normal.create(c(StDev = 99)))
# dist.normal<-dist.normal.create(c(0,1))
# summary(dist.normal)
# plot(dist.normal,plot=FALSE)
# plot(dist.normal,main="Normal",xlab="x",ylab="pdf")
# plot(dist.normal,main="Normal",what="CDF",xlab="x",ylab="cdf")
# plot(dist.normal,main="Normal",what="CDFinv",xlab="p",ylab="inv-cdf")
# plot(dist.normal,par=c(Mean=3))
# plot(dist.normal,par=c(StDev=3))

# dist.normal$par <- c(2,2)        # "permanent" Parameter parameter changenderung
# plot(dist.normal)
# plot(dist.normal,par=c(Mean=0))  # "temporary" parameter change

################################################################################

# Student distribution constructor:
# ---------------------------------
dist.student.create <- function(par=NA)
{
   # set default parameter values:

   par.default <- c( 0, 1, 3 )
   names(par.default) <- c( "Mean", "StDev", "DF" )

   # overwrite parameter values with provided values:

   p <- mergePar(par,par.default)

   # construct class:

   dist            <- list()
   dist$name       <- "Student"
   dist$range      <- function( par )
                      {
                         return( c( -Inf, +Inf ) )
                      }
   dist$par.names  <- names(p)
   dist$par.ranges <- matrix(c(-Inf, +Inf,  # Mean
                                  0, +Inf,   # Standard Deviation (StDev)
                                  3, +Inf),  # Degree of freedom (DF)
                      byrow=TRUE,ncol=2)
   dist$par        <- p
   dist$mean       <- function(par)
                      {
                         mean <- par[1]
                         names(mean) <- "Mean"
                         return(mean)
                       }
   dist$sd         <- function(par)
                      {
   	                  sd <- par[2]
   	                  names(sd) <- "StDev"
   	                  return(sd)
   	                  }
   dist$median     <- function(par)
                      {
                         median <- par[1]
                         names(median) <- "Median"
                         return(median)
                      }
   dist$mode       <- function(par)
                      {
                         mode <- par[1]
                         names(mode) <- "Mode"
                         return(mode)
                      }
   dist$pdf        <- function(x,par)
                      {
                         q       <- x
                         mu_t    <- par[1]
                         sigma_t <- par[2]
                         if(is.na(par[3])){N <- 3}
                         else {N <- par[3]}
                         return(dt((q-mu_t)/(sigma_t*sqrt((N-2)/N)), df = N)/(sigma_t*sqrt((N-2)/N)))
                      }
   dist$cdf        <- function(x,par)
                      {
                         q       <- x
                         mu_t    <- par[1]
                         sigma_t <- par[2]
                         if(is.na(par[3]))	{N <- 3}
                         else				{N <- par[3]}
                         return(pt( ((q-mu_t)/(sigma_t*sqrt((N-2)/N))),df = N))
                      }
   dist$cdf.inv    <- function(p,par)
                      {
                         mu_t    <- par[1]
                         sigma_t <- par[2]
                         if(is.na(par[3]))  {N <- 3}
                         else               {N <- par[3]}
                         return(qt(p,N)*sigma_t*sqrt((N-2)/N)+mu_t)
                      }
   class(dist)     <- "distribution"

   return(dist)
}

################################################################################

# Tests:
# ------
# print(dist.student.create(c(Mean=99)))
# print(dist.student.create(c(DF=99)))
# summary(dist.student.create(c(StDev=99)))
# dist.student <- dist.student.create(c(0,1,3))
# plot(dist.student,plot=FALSE)
# plot(dist.student,main="Student",xlab="x",ylab="pdf")
# plot(dist.student,main="Student",what="CDF",xlab="x",ylab="cdf")
# plot(dist.student,main="Student",what="CDFinv",xlab="p",ylab="inv-cdf")
# plot(dist.student,par=c(Mean=3))
# plot(dist.student,par=c(StDev=3))
# plot(dist.student,par=c(DF=3.5))


# dist.student$par <- c(2,2,10)    # "permanent" parameter change
# plot(dist.normal)
# plot(dist.normal,par=c(Mean=0))  # "temporary" parameter change

################################################################################

# Weibull distribution constructor:
# ---------------------------------

dist.weibull.create <- function(par=NA)
{
   # set default parameter values:

   par.default <- c(2,2)
   names(par.default) <- c( "Shape", "Scale" )

   # overwrite parameter values with provided values:

   p <- mergePar(par,par.default)

   # construct class:

   dist            <- list()
   dist$name       <- "Weibull"
   dist$range      <- function(par)
                      {
                         return(c(0,+Inf))
                      }
   dist$par.names  <- names(p)
   dist$par.ranges <- matrix(
                             c(1, +Inf,      # Shape
                               0, +Inf),     # Scale
                             byrow=TRUE,ncol=2)
   dist$par        <- p
   dist$mean       <- function(par)
                      {
                         mean <- par[2] * gamma(1+(1/par[1]))
   	                     names(mean) <- "Mean"
   	                     return(mean)
   	                  }
   dist$sd         <- function(par) 
                      {
                      	 sd <- ( par[2]^2 * gamma( 1+(2/par[1]) ) - (par[2] * gamma(1+(1/par[1])))^2 )^(1/2)  
                         names(sd) <- "StDev"
                         return(sd)
                      }
   dist$median     <- function(par)
                      { median <- par[2] *(log(2))^(1/par[1])
                      	names(median) <- "Median"
                      	return(median)
                      }
   dist$mode       <- function(par)
                      { 
                      	mode <- ifelse(par[1]>1, 
                                       par[2]*(((par[1]-1)/par[1])^(1/par[1])),
                                       NA)
                        names(mode) <- "Mode"
                      	return(mode) }
   dist$pdf        <- function(x,par) { return( dweibull(x, par[1], par[2]) ) }
   dist$cdf        <- function(x,par) { return( pweibull(x, par[1], par[2]) ) }
   dist$cdf.inv    <- function(p,par) { return( qweibull(p, par[1], par[2]) ) }
   class(dist)     <- "distribution"

   return(dist)
}

################################################################################

# Tests:
# ------
# print(dist.weibull.create(c(Shape=1.5)))
# print(dist.weibull.create(c(Scale=1.5)))
# dist.weibull <- dist.weibull.create(c(Shape=2,Scale=99))
# summary(dist.weibull)
# plot(dist.weibull,plot=FALSE)
# plot(dist.weibull,main="Weibull",xlab="x",ylab="pdf")
# plot(dist.weibull,main="Weibull",what="CDF",xlab="x",ylab="cdf")
# plot(dist.weibull,main="Weibull",what="CDFinv",xlab="p",ylab="inv-cdf")
# plot(dist.weibull,par=c(Shape=3))
# plot(dist.weibull,par=c(Scale=3))

# dist.weibull$par <- c(2,2)          # "permanent" parameter change
# plot(dist.weibull)
# plot(dist.weibull,par=c(Scale=20))  # "temporary" parameter change

################################################################################


# Lognormal distribution constructor:
# -----------------------------------

dist.lognormal.create <- function(par=NA)
{
   # set default parameter values:

   par.default <- c(1,1)
   names(par.default) <- c("Mean","StDev")

   # overwrite parameter values with provided values:

   p <- mergePar(par,par.default)

   # construct class:

   dist            <- list()
   dist$name       <- "Lognormal"
   dist$range      <- function(par)
                      {
                         return(c(0,+Inf))
                      }
   dist$par.names  <- names(p)
   dist$par.ranges <- matrix(
                             c(0, +Inf,      # mean
                               0, +Inf),     # sd
                             byrow=TRUE,ncol=2)
   dist$par        <- p
   dist$mean       <- function(par) { return(par[1]) }
   dist$sd         <- function(par) { return(par[2]) }
   dist$median     <- function(par) { mean   <- par[1]
   	                                  sd     <- par[2]
                                      median <- mean/(sqrt(1+(sd/mean)^2))
                                      names(median) <- "Median"
                                      return(median) }
   dist$mode       <- function(par) { mean <- par[1]
                                      sd   <- par[2]
                                      mode <- mean/(1+(sd/mean)^2)^(3/2)
                                      names(mode) <- "Mode"
                                      return(mode) }
   dist$pdf        <- function(x,par) { mean    <- par[1]
                                        sd      <- par[2]
                                        sdlog   <- sqrt(log(1+sd^2/mean^2))
                                        meanlog <- log(mean) - 0.5*sdlog^2
                                        return(dlnorm(x, meanlog, sdlog)) }
   dist$cdf        <- function(x,par) { mean    <- par[1]
                                        sd      <- par[2]
                                        sdlog   <- sqrt(log(1+sd^2/mean^2))
                                        meanlog <- log(mean) - 0.5*sdlog^2
                                        return(plnorm(x, meanlog, sdlog)) }
   dist$cdf.inv    <- function(p,par) { mean    <- par[1]
                                        sd      <- par[2]
                                        sdlog   <- sqrt(log(1+sd^2/mean^2))
                                        meanlog <- log(mean) - 0.5*sdlog^2
                                        return(qlnorm(p, meanlog, sdlog)) }
   class(dist)     <- "distribution"

   return(dist)
}

################################################################################

# Tests:
# ------
# print(dist.lognormal.create(c(2,1)))
# print(dist.lognormal.create(c(Mean=99)))
# print(dist.lognormal.create(c(StDev=99)))
# dist.lognormal <- dist.lognormal.create(c(1,1))
# summary(dist.lognormal)
# plot(dist.lognormal,plot=FALSE)
# plot(dist.lognormal,main="Lognormal",xlab="x",ylab="pdf")
# plot(dist.lognormal,main="Lognormal",what="CDF",xlab="x",ylab="cdf")
# plot(dist.lognormal,main="Lognormal",what="CDFinv",xlab="p",ylab="inv-cdf")
# plot(dist.lognormal,par=c(Mean=20))
# plot(dist.lognormal,par=c(Sd=1))
# dist.lognormal$par <- c(2,2)          # "permanent" parameter change
# plot(dist.lognormal)
# plot(dist.lognormal,par=c(StDev=20))  # "temporary" parameter change

################################################################################

# Beta distribution constructor:
# -----------------------------------

dist.beta.create <- function(par=NA)
{
   # set default parameter values:

   par.default <- c(1,1)
   names(par.default) <- c("Shape1","Shape2")

   # overwrite parameter values with provided values:

   p <- mergePar(par,par.default)

   # construct class:

   dist            <- list()
   dist$name       <- "Beta"
   dist$range      <- function(par)
                      {
                         return(c(0,1))
                      }
   dist$par.names  <- names(p)
   dist$par.ranges <- matrix(
                             c(1, +Inf,      # shape1
                               1, +Inf),     # shape2
                             byrow=TRUE,ncol=2)
   dist$par        <- p
   dist$mean       <- function(par) { mean <- par[1]/(par[1]+par[2]) 
   	                                  names(mean) <- "Mean"
   	                                  return(mean)}
   dist$sd         <- function(par) { sd <- sqrt(par[1]*par[2]/((par[1] + par[2])^2 * (par[1] + par[2] + 1)))
   	                                  names(sd) <- "StDev"
   	                                  return(sd) }
   dist$median     <- function(par) { median <- qbeta(0.5, par[1], par[2])
                                      names(median) <- "Median"
                                      return(median) }
   dist$mode       <- function(par) { 
                                      mode <- ifelse(par[1]>1&&par[2]>1,(par[1]-1)/(par[1]+par[2]-2),NA)
                                      names(mode) <- "Mode"
                                      return(mode) }
   dist$pdf        <- function(x,par) { 
                                        return(dbeta(x, par[1], par[2])) }
   dist$cdf        <- function(x,par) { 
                                        return(pbeta(x, par[1], par[2])) }
   dist$cdf.inv    <- function(p,par) { 
                                        return(qbeta(p, par[1], par[2])) }
   class(dist)     <- "distribution"

   return(dist)
}

################################################################################

# Tests:
# ------
# print(dist.beta.create(c(2,1)))
# print(dist.beta.create(c(shape1=99)))
# print(dist.beta.create(c(shape2=99)))
# dist.beta <- dist.beta.create(c(2,3))
# summary(dist.beta)
# plot(dist.beta,plot=FALSE)
# plot(dist.beta,main="beta",xlab="x",ylab="pdf")
# plot(dist.beta,main="beta",what="CDF",xlab="x",ylab="cdf")
# plot(dist.beta,main="beta",what="CDFinv",xlab="p",ylab="inv-cdf")
# plot(dist.beta,par=c(shape1=20))
# plot(dist.beta,par=c(shape2=20))
# dist.beta$par <- c(2,2)           # "permanent" parameter change
# plot(dist.beta)
# plot(dist.beta,par=c(shape1=20))  # "temporary" parameter change

################################################################################

# Gamma distribution constructor:
# -----------------------------------

dist.gamma.create <- function(par=NA)
{
   # set default parameter values:

   par.default <- c(1,1)
   names(par.default) <- c("shape","rate")

   # overwrite parameter values with provided values:

   p <- mergePar(par,par.default)

   # construct class:

   dist            <- list()
   dist$name       <- "Gamma"
   dist$range      <- function(par)
                      {
                         return(c(0,+Inf))
                      }
   dist$par.names  <- names(p)
   dist$par.ranges <- matrix(
                             c(0, +Inf,       # shape
                               0, +Inf),      # rate
                             byrow=TRUE,ncol=2)
   dist$par        <- p
   dist$mean       <- function(par) { mean <- par[1]*par[2] 
   	                                  names(mean) <- "Mean"
   	                                  return(mean)}
   dist$sd         <- function(par) { sd <- sqrt(par[1])*par[2]
   	                                  names(sd) <- "StDev"
   	                                  return(sd) }
   dist$median     <- function(par) { median <- qgamma(0.5, par[1], par[2])
                                      names(median) <- "Median"
                                      return(median) }
   dist$mode       <- function(par) { 
                                      mode <- ifelse(par[1]>=1,(par[1]-1)*(1/par[2]),NA)
                                      names(mode) <- "Mode"
                                      return(mode) }
   dist$pdf        <- function(x,par) { 
                                        return(dgamma(x, par[1], par[2])) }
   dist$cdf        <- function(x,par) { 
                                        return(pgamma(x, par[1], par[2])) }
   dist$cdf.inv    <- function(p,par) { 
                                        return(qgamma(p, par[1], par[2])) }
   class(dist)     <- "distribution"

   return(dist)
}

################################################################################

# Tests:
# ------
# print(dist.gamma.create(c(2,1)))
# print(dist.gamma.create(c(shape=99)))
# print(dist.gamma.create(c(rate=99)))
# dist.gamma <- dist.gamma.create(c(2,1))
# summary(dist.gamma)
# plot(dist.gamma,plot=FALSE)
# plot(dist.gamma,main="beta",xlab="x",ylab="pdf")
# plot(dist.gamma,main="beta",what="CDF",xlab="x",ylab="cdf")
# plot(dist.gamma,main="beta",what="CDFinv",xlab="p",ylab="inv-cdf")
# plot(dist.gamma,par=c(shape=20))
# plot(dist.gamma,par=c(rate=20))
# dist.gamma$par <- c(2,2)   # "permanent" parameter change
# plot(dist.gamma)
# plot(dist.gamma,par=c(shape=20))  # "temporary" parameter change

################################################################################


# F distribution constructor:
# -----------------------------------

dist.f.create <- function(par=NA)
{
   # set default parameter values:

   par.default <- c(3,5,0)
   names(par.default) <- c("df1","df2","ncp")

   # overwrite parameter values with provided values:

   p <- mergePar(par,par.default)

   # construct class:

   dist            <- list()
   dist$name       <- "F"
   dist$range      <- function(par)
                      {
                         return(c(0,+Inf))
                      }
   dist$par.names  <- names(p)
   dist$par.ranges <- matrix(
                             c(0, +Inf,       # df1
                               0, +Inf,       # df1
                            -Inf, +Inf),      # ncp
                             byrow=TRUE,ncol=2)
   dist$par        <- p
   dist$mean       <- function(par) { mean <- ifelse(par[2]>2,(par[2]*(par[1]+par[3]))/(par[1]*(par[2]-2)),NA)
   	                                  names(mean) <- "Mean"
   	                                  return(mean)}
   dist$sd         <- function(par) { sd <- ifelse(par[2]>4,
			sqrt(2*((par[2]/par[1])^2)*((((par[1]+par[3])^2)+(par[1]+2*par[3])*(par[2]-2))/((((par[2]-2)^2)*(par[2]-4))))),
			NA)
   	                                  names(sd) <- "StDev"
   	                                  return(sd) }
   dist$median     <- function(par) { median <- qf(0.5, par[1], par[2], par[3])
                                      names(median) <- "Median"
                                      return(median) }
   dist$mode       <- function(par) { #m<-which.max(plot(dist.f,par=c(df1=par[1],df2=par[2],ncp=par[3]),plot=FALSE)[,2])
                                      #mode <-      plot(dist.f,par=c(df1=par[1],df2=par[2],ncp=par[3]),plot=FALSE)[m,1]
                                      #names(mode) <- "Approximated Mode"
                                      if(par[2]>2){mode <- ((par[1]-2)/par[1])*(par[2]/(par[2]+2))}
                                      if(par[2]<=2){mode <- NA}
                                      names(mode) <- "Mode"
                                      return(mode) }
   dist$pdf        <- function(x,par) { 
                                        return(df(x, par[1], par[2], par[3])) }
   dist$cdf        <- function(x,par) { 
                                        return(pf(x, par[1], par[2], par[3])) }
   dist$cdf.inv    <- function(p,par) { 
                                        return(qf(p, par[1], par[2], par[3])) }
   class(dist)     <- "distribution"

   return(dist)
}

################################################################################

# Tests:
# ------
# print(dist.f.create(c(2,3,0)))
# print(dist.f.create(c(df1=99)))
# print(dist.f.create(c(df2=99)))
# print(dist.f.create(c(ncp=99)))
# dist.f <- dist.f.create(c(3,5,0))
# summary(dist.f)
# plot(dist.f,plot=FALSE)
# plot(dist.f,main="beta",xlab="x",ylab="pdf")
# plot(dist.f,main="beta",what="CDF",xlab="x",ylab="cdf")
# plot(dist.f,main="beta",what="CDFinv",xlab="p",ylab="inv-cdf")
# plot(dist.f,par=c(df2=20))
# plot(dist.f,par=c(df1=2))
# dist.f$par <- c(2,2,0)   # "permanent" parameter change
# plot(dist.f)
# plot(dist.f,par=c(df2=20))  # "temporary" parameter change

################################################################################

# Uniform distribution constructor:
# -----------------------------------

dist.uniform.create <- function(par=NA)
{
   # set default parameter values:

   par.default <- c(0,1)
   names(par.default) <- c("Min","Max")

   # overwrite parameter values with provided values:

   p <- mergePar(par,par.default)

   # construct class:

   dist            <- list()
   dist$name       <- "Uniform"
   dist$range      <- function(par)
                      {
                         return(c(par[1],par[2]))
                      }
   dist$par.names  <- names(p)
   dist$par.ranges <- matrix(
                             c(-Inf, +Inf,       # Min
                               -Inf, +Inf),      # Max
                             byrow=TRUE,ncol=2)
   dist$par        <- p
   dist$mean       <- function(par) { mean <- (par[1]+par[2])/2
   	                                  names(mean) <- "Mean"
   	                                  return(mean)}
   dist$sd         <- function(par) { sd <- sqrt((par[2]-par[1])^2/12)
   	                                  names(sd) <- "StDev"
   	                                  return(sd) }
   dist$median     <- function(par) { median <- (par[1]+par[2])/2
                                      names(median) <- "Median"
                                      return(median) }
   dist$mode       <- function(par) { 
                                      mode <- NA 
                                      names(mode) <- "Mode"
                                      return(mode) }
   dist$pdf        <- function(x,par) { 
                                        return(dunif(x, par[1], par[2])) }
   dist$cdf        <- function(x,par) { 
                                        return(punif(x, par[1], par[2])) }
   dist$cdf.inv    <- function(p,par) { 
                                        return(qunif(p, par[1], par[2])) }
   class(dist)     <- "distribution"

   return(dist)
}

################################################################################

# Tests:
# ------
# print(dist.uniform.create(c(1,99)))
# print(dist.uniform.create(c(Min=99)))
# print(dist.uniform.create(c(Max=99)))
# dist.uniform <- dist.uniform.create(c(-1,1))
# summary(dist.uniform)
# plot(dist.uniform,plot=FALSE)
# plot(dist.uniform,main="uniform",xlab="x",ylab="pdf")
# plot(dist.uniform,main="beta",what="CDF",xlab="x",ylab="cdf")
# plot(dist.uniform,main="beta",what="CDFinv",xlab="p",ylab="inv-cdf")
# plot(dist.uniform,par=c(Max=20))
# plot(dist.uniform,par=c(Min=2))
# dist.uniform$par <- c(2,3)   # "permanent" parameter change
# plot(dist.uniform)
# plot(dist.uniform,par=c(Max=20))  # "temporary" parameter change


################################################################################

# Logistic distribution constructor:
# -----------------------------------

dist.logistic.create <- function(par=NA)
{
   # set default parameter values:

   par.default <- c(0,1)
   names(par.default) <- c("Location","Scale")

   # overwrite parameter values with provided values:

   p <- mergePar(par,par.default)

   # construct class:

   dist            <- list()
   dist$name       <- "Logistic"
   dist$range      <- function(par)
                      {
                         return(c(-Inf,+Inf))
                      }
   dist$par.names  <- names(p)
   dist$par.ranges <- matrix(
                             c(-Inf, +Inf,       # Location
                                  0, +Inf),      # Scale
                             byrow=TRUE,ncol=2)
   dist$par        <- p
   dist$mean       <- function(par) { mean <- par[1]
   	                                  names(mean) <- "Mean"
   	                                  return(mean)}
   dist$sd         <- function(par) { sd <- sqrt((pi^2/3)*par[2]^2)
   	                                  names(sd) <- "StDev"
   	                                  return(sd) }
   dist$median     <- function(par) { median <- par[1]
                                      names(median) <- "Median"
                                      return(median) }
   dist$mode       <- function(par) { 
                                      mode <- par[1] 
                                      names(mode) <- "Mode"
                                      return(mode) }
   dist$pdf        <- function(x,par) { return(dlogis(x, par[1], par[2])) }
   dist$cdf        <- function(x,par) { return(plogis(x, par[1], par[2])) }
   dist$cdf.inv    <- function(p,par) { return(qlogis(p, par[1], par[2])) }
   class(dist)     <- "distribution"

   return(dist)
}

################################################################################

# Tests:
# ------
# print(dist.logistic.create(c(2,1)))
# print(dist.logistic.create(c(Location=99)))
# print(dist.logistic.create(c(Scale=99)))
# dist.logistic <- dist.logistic.create(c(1,3))
# summary(dist.logistic)
# plot(dist.logistic,plot=FALSE)
# plot(dist.logistic,main="beta",xlab="x",ylab="pdf")
# plot(dist.logistic,main="beta",what="CDF",xlab="x",ylab="cdf")
# plot(dist.logistic,main="beta",what="CDFinv",xlab="p",ylab="inv-cdf")
# plot(dist.logistic,par=c(Location=20))
# plot(dist.logistic,par=c(Scale=2))
# dist.logistic$par <- c(2,3)   # "permanent" parameter change
# plot(dist.logistic)
# plot(dist.logistic,par=c(Location=20))  # "temporary" parameter change

################################################################################

# Required elements of class "transformations":
# ------------------------------------------
#
# name            character        name
# range.x         function         range over which the transformation is defined
# range.y         function         range of the domain of the transformations image
# par.names       character(npar)  parameter names (npar names)
# par.ranges      matrix(npar,2)   ranges of parameters (one after the other)
# par             numeric(npar)    (default) parameter values (npar values)
# trans.forw      function         forward transformation
# trans.backw     function         backward transformation
# trans.deriv     function         derivation of transformation


# Methods of class "transformation":
# --------------------------------
#
# print(trans,...)
# summary(trans,...)
# plot(trans, par = trans$par, range = NA, plot = TRUE,...)     produces a pdf listing if plot=FALSE

print.transformation <- function(x = trans,...)
{  trans <- x
   rm(x)
   title <- paste(trans$name,"Transformation")
   cat(paste(title,"\n"))
   cat(paste(paste(rep("-",nchar(title)),collapse=""),"\n"))
   cat(paste("Range of x:",trans$range.x(trans$par)[1],trans$range.x(trans$par)[2],"\n"))
   cat(paste("Range of y:",trans$range.y(trans$par)[1],trans$range.y(trans$par)[2],"\n"))
   cat("Parameter names, ranges and values:\n")
   ranges <- trans$par.ranges
   tab <- cbind(ranges,trans$par)
   colnames(tab) <- c("Min","Max","Value")
   rownames(tab) <- trans$par.names
   print(tab)
   cat("\n")
}

summary.transformation <- function(object,...)
{
   print(object,...)
}

plot.transformation <- function(x = trans, par = trans$par, range.x = NA, range.y = NA, what = "TRANS.FORW", plot = TRUE, length = 101,...)
{  trans <- x
   rm(x)
   # theoretic limits
   x.min <- as.numeric(trans$range.x(par)[1])
   x.max <- as.numeric(trans$range.x(par)[2])
   y.min <- as.numeric(trans$range.y(par)[1])
   y.max <- as.numeric(trans$range.y(par)[2])
   if ( what == "TRANS.FORW" )
   {
      if ( is.na( range.x[1] ) && is.na( range.y[1]) )
      {  
         if( x.min>-Inf && x.max<Inf )
         {  
            diff.x <- 0.01*(x.max-x.min)
            x.min  <- x.min+diff.x
            x.max  <- x.max-diff.x
         }
         if( y.min>-Inf && y.max<Inf )
         {  
            diff.y <- 0.01*(y.max-y.min)
            x.min  <- TRANS.BACKW(trans,y.min+diff.y,par)
            x.max  <- TRANS.BACKW(trans,y.max-diff.y,par)
         }
         if(x.min>-Inf && x.max<Inf && y.min>-Inf && y.max<Inf )
         {
            diff.x <- 0.01*(x.max-x.min)
            diff.y <- 0.01*(y.max-y.min)
            x.min <- min(x.min+diff.x,TRANS.BACKW(trans,y.min+diff.y,par))
            x.max <- max(x.max-diff.x,TRANS.BACKW(trans,y.max-diff.y,par))
         }
         if( !  ((x.min>-Inf && x.max<Inf )|(y.min>-Inf && y.max<Inf )) )
         {
            stop("*** 'range.x  = c( ., .)' must be specifyed in the function plot.transformation TRANS.FORW")
         }
      }
      else
      {
      	 if( (!(is.na(range.x[1]) && is.na(range.x[2]))) && (is.na(range.y[1]) && is.na(range.y[2])) )
      	 {
            x.min <- range.x[1]
            x.max <- range.x[2]
         }
         if(! (is.na(range.y[1]) && is.na(range.y[2])) && (is.na(range.x[1]) && is.na(range.x[2])))
         {
            x.min <- TRANS.BACKW(trans,range.y[1],par)
            x.max <- TRANS.BACKW(trans,range.y[2],par)
         }
         if((!(is.na(range.x[1]) && is.na(range.x[2]))) && (!(is.na(range.y[1]) && is.na(range.y[2]))) )
         {
            x.min <- min(range.x[1],TRANS.BACKW(trans,range.y[1],par))
            x.max <- max(range.x[2],TRANS.BACKW(trans,range.y[2],par))
         }
      }
      x <- seq(x.min,x.max,length=length)
      y <- TRANS.FORW(trans,x,par)
   }
   if ( what == "TRANS.BACKW" )
   {
      if ( is.na( range.x[1] ) && is.na( range.y[1]) )
      {  
         if( x.min>-Inf && x.max<Inf )
         {  
            diff.x <- 0.01*(x.max-x.min)
            y.min  <- TRANS.FORW(trans,x.min+diff.x,par)
            y.max  <- TRANS.FORW(trans,x.max-diff.x,par)
         }
         if( y.min>-Inf && y.max<Inf )
         {  
            diff.y <- 0.01*(y.max-y.min)
            y.min  <- y.min+diff.y
            y.max  <- y.max-diff.y
         }
         if(x.min>-Inf && x.max<Inf && y.min>-Inf && y.max<Inf )
         {
            diff.x <- 0.01*(x.max-x.min)
            diff.y <- 0.01*(y.max-y.min)
            y.min <- min(TRANS.FORW(trans,x.min+diff.x,par),y.min+diff.y)
            y.max <- max(TRANS.FORW(trans,x.max-diff.x,par),y.max-diff.y)
         }
         if( !  ((x.min>-Inf && x.max<Inf )|(y.min>-Inf && y.max<Inf )) )
         {
            stop("*** 'range.y = c( ., .)' must be specifyed in the function plot.transformation TRANS.BACKW")
         }
      }
      else
      {
      	 if( (!(is.na(range.x[1]) && is.na(range.x[2]))) && (is.na(range.y[1]) && is.na(range.y[2])))
      	 {
            y.min <- TRANS.FORW(trans,range.x[1],par)
            y.max <- TRANS.FORW(trans,range.x[2],par)
         }
         if(! (is.na(range.y[1]) && is.na(range.y[2])) && (is.na(range.x[1]) && is.na(range.x[2])))
         {
            y.min <- range.y[1]
            y.max <- range.y[2]
         }
         if((!(is.na(range.x[1]) && is.na(range.x[2]))) && (!(is.na(range.y[1]) && is.na(range.y[2]))) )
         {
            y.min <- min(TRANS.FORW(trans,range.x[1],par),range.y[1])
            y.max <- max(TRANS.FORW(trans,range.x[2],par),range.y[2])
         }
      }
      x <- seq(y.min,y.max,length=length)
      y <- TRANS.BACKW(trans,x,par)
   }
   if ( what == "TRANS.DERIV" )
   {
      if ( is.na( range.x[1] ) && is.na( range.y[1]) )
      {  
         if( x.min>-Inf && x.max<Inf )
         {  
            diff.x <- 0.05*(x.max-x.min)
            x.min  <- x.min+diff.x
            x.max  <- x.max-diff.x
         }
         if( y.min>-Inf && y.max<Inf )
         {  
            diff.y <- 0.05*(y.max-y.min)
            x.min  <- TRANS.BACKW(trans,y.min+diff.y,par)
            x.max  <- TRANS.BACKW(trans,y.max-diff.y,par)
         }
         if(x.min>-Inf && x.max<Inf && y.min>-Inf && y.max<Inf )
         {
            diff.x <- 0.05*(x.max-x.min)
            diff.y <- 0.05*(y.max-y.min)
            x.min <- min(x.min+diff.x,TRANS.BACKW(trans,y.min+diff.y,par))
            x.max <- max(x.max-diff.x,TRANS.BACKW(trans,y.max-diff.y,par))
         }
         if( !  ((x.min>-Inf && x.max<Inf )|(y.min>-Inf && y.max<Inf )) )
         {
            stop("*** 'range.x = c(.,.)' must be specifyed in the function plot.transformation TRANS.DERIV")
         }
      }
      else
      {
      	 if( (!(is.na(range.x[1]) && is.na(range.x[2]))) && (is.na(range.y[1]) && is.na(range.y[2])) )
      	 {
            x.min <- range.x[1]
            x.max <- range.x[2]
         }
         if(! (is.na(range.y[1]) && is.na(range.y[2])) && (is.na(range.x[1]) && is.na(range.x[2])) )
         {
            stop("*** 'range.x = c(.,.)' must be specifyed in the function plot.transformation TRANS.DERIV")
         }
         if( (!(is.na(range.x[1]) && is.na(range.x[2]))) && (!(is.na(range.y[1]) && is.na(range.y[2]))) )
         {
            x.min <- range.x[1]
            x.max <- range.x[2]
         }
      }
      x <- seq(x.min,x.max,length=length)
      y <- TRANS.DERIV(trans,x,par)
   }
   if ( !what == "TRANS.FORW"  && !what == "TRANS.BACKW"  && !what == "TRANS.DERIV" )
   {
         stop("*** illegal value of 'what' in plot.transformation")
   }   
   if ( plot )
   {
   	  res <- cbind(x,y)
   	  x <- res[,1]
   	  y <- res[,2]
      plot(x, y, type="l" , ...)
   }
   if ( ! plot )
   {
      res <- cbind( x, y)
      colnames(res) <- c( "x", "y" )
      return(res)
   }
}


trans.arctan.create <- function(par=NA)
{
   # set default parameter values:

   par.default <- c(0,1)
   names(par.default) <- c("Min","Max")

   # overwrite parameter values with provided values:

   p <- mergePar(par,par.default)

   # construct class:

   trans             <- list()
   trans$name        <- "Arctan"
   trans$range.x     <- function(par){return(c(-Inf,+Inf))}
   trans$range.y     <- function(par){return(par)}
   trans$par.names   <- names(p)
   trans$par.ranges  <- matrix(
                               c(-Inf, +Inf,      # Min
                                 -Inf, +Inf),     # Max
                               byrow=TRUE,ncol=2)
   trans$par         <- p
   trans$trans.forw  <- function(x,par)
                        { y <- 0.5*(par[1]+par[2]) + (par[2]-par[1])/pi*atan(x)
                          return(as.numeric(y)) }
   trans$trans.backw <- function(y,par)
                        { x <- tan(0.5*pi*(2*y-par[2]-par[1])/(par[2]-par[1]))
                          return(as.numeric(x)) }
   trans$trans.deriv <- function(x,par)
                        { dydx <- (par[2]-par[1])/pi/(1+x*x)
                          return(as.numeric(dydx)) }
   class(trans)      <- "transformation"

   return(trans)
}
# Tests
# trans.arctan <- trans.arctan.create(c(0,10))
# x11()
# plot(trans.arctan)
# plot(trans.arctan,what = "TRANS.BACKW")
# plot(trans.arctan,what = "TRANS.DERIV")


trans.tan.create <- function(par=NA)
{
   # set default parameter values:

   par.default <- c(0,1)
   names(par.default) <- c("Min","Max")

   # overwrite parameter values with provided values:

   p <- mergePar(par,par.default)

   # construct class:

   trans             <- list()
   trans$name        <- "Tan"
   trans$range.x     <- function(par){return(par)}
   trans$range.y     <- function(par){return(c(-Inf,Inf))}
   trans$par.names   <- names(p)
   trans$par.ranges  <- matrix(
                               c(-Inf, +Inf,      # Min
                                 -Inf, +Inf),     # Max
                               byrow=TRUE,ncol=2)
   trans$par         <- p
   trans$trans.forw  <- function(x,par)
                        { y <- tan(0.5*pi*(2*x-par[2]-par[1])/(par[2]-par[1]))
                          return(as.numeric(y)) }
   trans$trans.backw <- function(y,par)
                        { x <- 0.5*(par[1]+par[2]) + (par[2]-par[1])/pi*atan(y)
                          return(as.numeric(x)) }
   trans$trans.deriv <- function(x,par)
                        { dydx <- pi/(par[2]-par[1])/
                                  (cos(0.5*pi*(2*x-par[2]-par[1])/
                                       (par[2]-par[1]))^2)
                          return(as.numeric(dydx)) }
   class(trans)      <- "transformation"

   return(trans)
}
# Tests
# trans.tan <- trans.tan.create(c(0,10))
# x11()
# plot(trans.tan)
# plot(trans.tan,what = "TRANS.BACKW")
# plot(trans.tan,what = "TRANS.DERIV")


trans.log.create <- function(par=NA)
{
   # set default parameter values:

   par.default <- c(NA)
   names(par.default) <- c("-")

   # overwrite parameter values with provided values:

   p <- mergePar(par,par.default)

   # construct class:

   trans             <- list()
   trans$name        <- "Log"
   trans$range.x     <- function(par){return(c(0.000001,Inf))}
   trans$range.y     <- function(par){return(c(-Inf,Inf))}
   trans$par.names   <- names(p)
   trans$par.ranges  <- matrix(c(NA),     
                               byrow=TRUE,ncol=2)
   trans$par         <- p
   trans$trans.forw  <- function(x,par)
                        { y <- log(x)
                          return(as.numeric(y)) }
   trans$trans.backw <- function(y,par)
                        { x <- exp(y)
                          return(as.numeric(x)) }
   trans$trans.deriv <- function(x,par)
                        { dydx <- 1/x
                          return(as.numeric(dydx)) }
   class(trans)      <- "transformation"

   return(trans)
}
# Tests
# trans.log <- trans.log.create()
# x11()
# plot(trans.log,range.x=c(-1,1))
# plot(trans.log,what = "TRANS.BACKW",range.y=c(-1,1))
# plot(trans.log,what = "TRANS.DERIV",range.x=c(-1,1))


trans.dil.create <- function(par=NA)
{
   # set default parameter values:

   par.default <- c(0,1,0,1)
   names(par.default) <- c("Min1","Max1","Min2","Max2")

   # overwrite parameter values with provided values:

   p <- mergePar(par,par.default)

   # construct class:

   trans             <- list()
   trans$name        <- "Dilation"
   trans$range.x     <- function(par){return(c(par[1],par[2]))}
   trans$range.y     <- function(par){return(c(par[3],par[4]))}
   trans$par.names   <- names(p)
   trans$par.ranges  <- matrix(
                               c(-Inf, +Inf,      # Min1
                                 -Inf, +Inf,      # Max1
                                 -Inf, +Inf,      # Min2
                                 -Inf, +Inf),     # Max2
                               byrow=TRUE,ncol=2)
   trans$par         <- p
   trans$trans.forw  <- function(x,par)
                        { y <- (x-par[1]) * (par[4]-par[3])/(par[2]-par[1]) + par[3]
                          return(as.numeric(y)) }
   trans$trans.backw <- function(y,par)
                        { x <- (y-par[3]) * (par[2]-par[1])/(par[4]-par[3]) + par[1]
                          return(as.numeric(x)) }
   trans$trans.deriv <- function(x,par)
                        { dydx <- (par[4]-par[3])/(par[2]-par[1])
                          return(as.numeric(dydx)) }
   class(trans)      <- "transformation"
   return(trans)
}
# Tests
# trans.dil <- trans.dil.create(c(0,1,4,5))
# x11()
# plot(trans.dil,range.x=c(-1,1))
# plot(trans.dil,what = "TRANS.BACKW", range.y = c(-1,1))
# plot(trans.dil,what = "TRANS.DERIV", range.x = c(-1,1))


trans.exp.create <- function(par=NA)
{
   # set default parameter values:

   par.default <- c(0,1,0)
   names(par.default) <- c("a","b","c")

   # overwrite parameter values with provided values:

   p <- mergePar(par,par.default)

   # construct class: 

   trans             <- list()
   trans$name        <- "Exponential"
   trans$range.x     <- function(par){return( c( -Inf, Inf ))}
   trans$range.y     <- function(par){return( c( -Inf, Inf ))}
   trans$par.names   <- names(p)
   trans$par.ranges  <- matrix(
                               c( 0, +Inf,      # a
                                  0, +Inf,      # b
                                  0, +Inf),     # c
                               byrow=TRUE,ncol=2)
   trans$par         <- p
   trans$trans.forw  <- function( x, par )
                        { y <- -(par[1]/par[2]^2)*exp(-par[2]*x)+par[3]*x+(par[1]/par[2]^2)
                          return(as.numeric(y)) }
   trans$trans.backw <- function( y, par )
                        {   
                        	a <- as.numeric(par[1])
                        	b <- as.numeric(par[2])
                        	c <- as.numeric(par[3])
                        	x <- y
                        	for(i in 1:length(y))
                        	{
		                       x[i] <- uniroot(function(x)
		                       {
		                          -(a/b^2)*exp(-b*x) + c*x + (a/b^2) - y[i]
		                       	},
		                       lower   = -999999999999999999999999,
		                       upper   =  999999999999999999999999,
		                       f.lower = -999999999999999999999999,
		                       f.upper =  999999999999999999999999)$root
		                    }
                          return(as.numeric(x)) }
   trans$trans.deriv <- function( x, par )
                        { dydx <- (par[1]/par[2])*exp(-par[2]*x)+par[3]
                          return(as.numeric(dydx)) }
   class(  trans )      <- "transformation"
   return( trans )
}

# Tests
# trans.exp <- trans.exp.create(c(3, 2, 1))
# x11()
# plot(trans.exp,range.x=c(-1,1))
# plot(trans.exp,what = "TRANS.BACKW", range.y = c(-4,3))
# plot(trans.exp,what = "TRANS.DERIV", range.x = c(-1,1))


TRANS.RANGE.X.transformation <- function(trans, par = trans$par,...)
{
   par.local <- mergePar(par,trans$par,trans$par.names)
   return(trans$range.x(par.local))
}

TRANS.RANGE.Y.transformation <- function(trans, par = trans$par,...)
{
   par.local <- mergePar(par,trans$par,trans$par.names)
   return(trans$range.y(par.local))
}

TRANS.FORW.transformation <- function(trans,x, par = trans$par,...)
{
   par.local <- mergePar(par,trans$par,trans$par.names)
   return(trans$trans.forw(x,par.local))
}

TRANS.BACKW.transformation <- function(trans,y, par = trans$par,...)
{
   par.local <- mergePar(par,trans$par,trans$par.names)
   return(trans$trans.backw(y,par.local))
}

TRANS.DERIV.transformation <- function(trans, x, par = trans$par,...)
{
   par.local <- mergePar(par,trans$par,trans$par.names)
   return(trans$trans.deriv(x,par.local))
}

TRANS.RANGE.X <- function(trans, par = trans$par,...) UseMethod("TRANS.RANGE.X")
TRANS.RANGE.Y <- function(trans, par = trans$par,...) UseMethod("TRANS.RANGE.Y")
TRANS.FORW  <- function(trans,x, par = trans$par,...) UseMethod("TRANS.FORW")
TRANS.BACKW <- function(trans,y, par = trans$par,...) UseMethod("TRANS.BACKW")
TRANS.DERIV <- function(trans,x, par = trans$par,...) UseMethod("TRANS.DERIV")

# creator for a transformed distribution

dist.trans.create <- function(dist,trans)
{
   # check consistency:
   
   if ( dist$range(dist$par)[1] != trans$range.y(trans$par)[1] )
   {
      stop("*** lower range of the variable to be tranformed not equal to lower range of distribution")
   }
   if ( dist$range(dist$par)[2] != trans$range.y(trans$par)[2] )
   {
      stop("*** upper range of the variable to be tranformed not equal to lower range of distribution")
   }
   
   # construct class:
   
   dist.trans            <- list()
   dist.trans$dist       <- dist
   dist.trans$trans      <- trans
   dist.trans$name       <- paste(dist$name,trans$name,sep="_")
   dist.trans$range      <- function(par){ return( trans$range.y( par ) ) }
   dist.trans$par.names  <- c(dist$par.names,trans$par.names)
   dist.trans$par.ranges <- rbind(dist$par.ranges,trans$par.ranges)
   dist.trans$par        <- c(dist$par,trans$par)
   dist.trans$mean       <- function(par=NA) { return(NA) }
   dist.trans$sd         <- function(par=NA) { return(NA) }
   dist.trans$median     <- function(par=NA) { return(NA) }
   dist.trans$mode       <- function(par=NA) { return(NA) }
   dist.trans$pdf        <- function(x,par)  { return(PDF(dist,x,par)) }
   dist.trans$cdf        <- function(x,par)  { return(CDF(dist,x,par)) } 
   dist.trans$cdf.inv    <- function(p,par)  { return(CDFinv(dist,p,par)) }
   class(dist.trans)     <- "distribution"
   
   return(dist.trans)
}

################################################################################

# Elements of class "drclass":
# ---------------------------------
#
# name            character        name
# range           vector           returns the range over the drc is defined
# p               vector           probabilities
# q               vector           quantiles
# dist.lower      object           object in the class "distribution"
# dist.upper      object           object in the class "distribution"
# alpha           numeric          confidence level for metrics

################################################################################

# DRC constructor:
# ----------------

drclass.create <- function(p = c(0.05,0.25,0.5,0.75,0.95),
                           q = qnorm(c(0.05,0.25,0.5,0.75,0.95)),
                           dist.lower = dist.normal.create(c(0,1)),
                           dist.upper = dist.normal.create(c(0,1))  )
{

drc                      <- list()
drc$name                 <- "Default"
drc$range                <- dist.upper$range(dist.upper$par) 
drc$p                    <- p
drc$q                    <- q
drc$dist.lower           <- dist.lower
drc$dist.upper           <- dist.upper
class(drc)               <- "drclass"
return(drc)

}

# Methods of class "drclass":
# --------------------------------
#
# print(drc, ...)                       prints the drc without the metrics (see summary).
# summary(drc, alpha, ...)              prints the drc attributes and the metrics for alpha.
# plot(drc,range=NA, , plot.stat.values = FALSE, makePDF = FALSE, ...)
#                                       plot.stat.values: plots the means, sds, medians and modes.
#                                       makePDF: if TRUE creates a sketch.pdf
# metric.ci(drc, alpha, ...)            returns a vector of the credib. interval / conf. int.
# metric.width(drc, alpha, ...)         returns the metric of width
# metric.shape(drc, alpha, ...)         returns the metric of shape
# metric.mode(drc, alpha, ...)          returns the metric of mode
# metrics(drc, alpha, ...)              returns the confidence interval and the 3 metrics
# Fl.drclass(drc, x, ...)               returns lower envelope of CDFs
# FlInv(drc, x, ...)                    returns inv. lower envelope of CDFs
# Fu(drc, x, ...)                       returns upper envelope of CDFs
# FuInv(drc, x, ...)                    returns inv. upper envelope of CDFs
# Kappa(drc,...)                        calculates and returns Kappa
# Lambda(drc,...)                       calculates and returns Lambda
print.drclass <- function(x = drc,...)
{  drc <- x
   rm(x)
   title <- paste(drc$name," Density Ratio Class")
   cat(paste(paste(rep("-",nchar(title)),collapse=""),"\n"))
   cat(paste(title,"\n"))
   cat(paste(paste(rep("-",nchar(title)),collapse=""),"\n"))
   cat(paste("\n ----------------- \n"))
   cat(paste("\n ----- Range ----- \n"))
   cat(paste("\n ----------------- \n"))
   values <- matrix(c(  drc$range[1],
   						drc$range[2]
   						),ncol=2,byrow=TRUE)
   colnames(values) <- c("Min","Max")
   rownames(values) <- c("Range")
   print(values)
   cat(paste("\n ----------------- \n"))
   cat(paste("\n ----- P & Q ----- \n"))
   cat(paste("\n ----------------- \n"))
   values <- matrix(c(  drc$p,
   						drc$q
   						),ncol=2,byrow=FALSE)
   rownames(values) <- rep("Pt.",length(values[,1]))
   colnames(values) <- c("Probabilities",
   						 "Quantiles")
   print(values)  
   # cat(paste("\n ----- Transformation ----- \n"))
   # cat("...none implemented so far...\n")
   cat(paste("\n ----------------- \n"))
   cat(paste("\n ----- Lower ----- \n"))
   cat(paste("\n ----------------- \n"))
   cat(paste(print(drc$dist.lower),"\n"))
   cat(paste("\n ----------------- \n"))
   cat(paste("\n ----- Upper ----- \n"))
   cat(paste("\n ----------------- \n"))
   cat(paste(print(drc$dist.upper),"\n"))
   cat(paste("\n -------------------------- \n"))
   cat(paste("\n ----- Kappa & Lambda ----- \n"))
   cat(paste("\n -------------------------- \n"))
   values <- matrix(c(  Lambda(drc),
   						Kappa(drc)
   						),ncol=1,byrow=TRUE)
   colnames(values) <- c("Value")
   rownames(values) <- c("Lambda",
   						 "Kappa")
   print(values)
   cat("\n")
}

summary.drclass <- function(object = drc, alpha = 0.05, ...)
{  drc <- object
   rm(object)
   print(drc)
   cat(paste("\n ------------------- \n"))
   cat(paste("\n ----- Metrics ----- \n"))
   cat(paste("\n ------------------- \n"))
   cat(paste("Confidence level: alpha = ",alpha,"\n"))
   cat(paste("Confidence interval bound left:  ",metric.ci(drc,alpha)[1],"\n"))
   cat(paste("Confidence interval bound right: ",metric.ci(drc,alpha)[2],"\n"))
   cat(paste("Metric of width:                 ",metric.width(drc,alpha),"\n"))
   cat(paste("Metric of shape:                 ",metric.shape(drc,alpha),"\n"))
   cat(paste("Metric of mode:                  ",metric.mode(drc,alpha),"\n"))
   cat("\n")
}

plot.drclass <- function(x = drc, range = NA, plot.stat.values = FALSE, makePDF = FALSE, ...)
{   drc <- x
	rm(x)
	p <- drc$p
	q <- drc$q
	k <- Kappa(drc)
	
	if( is.na(range)[1] )
	{
		if( drc$range[1] > -Inf && drc$range[2] < Inf )
		{
			r <- drc$range
		}
		else
        {
		    r <- c(FuInv(drc,0.01), FlInv(drc,0.99))
		}
	}
	
	if( ! is.na(range)[1] )
	{
		r <- range
	}
	
	if(makePDF)
	{
	   pdf.height 	<-  9
	   pdf.width  	<- 16
	   # Mac specific with date and time:
	   # pdf(paste("sketch_",round(k,2),"_",format(Sys.time(), "%a %b %d %H:%M:%S %Y"),".pdf",sep =""),
	   pdf(paste("sketch_",round(k,2),".pdf",sep =""),
       height=pdf.height,width=pdf.width)
    }
    
	par(mfrow=c(1,2))
	
	plot(numeric(0),numeric(0),mar=c(5, 4, 4, 11)+.1, xpd=TRUE,
		title("DRC Envelope of all CDFs",
		sub = paste("lower dist",drc$dist.lower$name,
				" upper dist",drc$dist.upper$name,
				" with Kappa =",round(k,2)),
    	cex.main = 1.5, font.main= 2, col.main= "black",
    	cex.sub = 1, font.sub = 3, col.sub = "black"),
		xlim = r,
		ylim = c(0,1),
		xlab = (expression(theta)),
		ylab ="Probability",
		cex.main = 1.5,
		cex.axis = 1.5,
		cex.lab  = 1.5)

		points(q, p, lty = 2, lwd = 1)
		for(i in 1:length(p)){
			#lines(c(FuInv(p[i],qmin,qmax,l,u,k,par_l,par_u,par2)$root,qmin[i]),
			#	  c(p[i],p[i]),lty=3,lwd=4)
			#lines(c(FlInv(p[i],qmin,qmax,l,u,k,par_l,par_u,par2)$root,qmax[i]),
			#	  c(p[i],p[i]),lty=3,lwd=4) 
		}
		
	    curve(Fl(drc,x),
		   min(r),
		   max(r),
		   col=1,
		   add = TRUE,
		   lty="solid",
		   lwd=2)
	    curve(Fu(drc,x),
		   min(r),
		   max(r),
		   col=1,
		   add = TRUE,
		   lty="solid",
		   lwd=2)

	xlpdf <- plot(drc$dist.lower, range = r, what="PDF", plot=FALSE, length=101)[,1]
	ylpdf <- plot(drc$dist.lower, range = r, what="PDF", plot=FALSE, length=101)[,2]
	
	xupdf <- plot(drc$dist.upper, range = r, what="PDF", plot=FALSE, length=101)[,1]
	yupdf <- plot(drc$dist.upper, range = r, what="PDF", plot=FALSE, length=101)[,2]
    
	# DRC PDFs
	plot(numeric(0),numeric(0),
		title("DRC PDFs",
		sub = paste("lower dist",drc$dist.lower$name,"  upper dist",drc$dist.upper$name," with Kappa =",round(k,2)),
		cex.main = 1.5,   font.main = 2,    col.main = "black",
		cex.sub  = 1,     font.sub  = 3,    col.sub  = "black"),
		xlim = r,
		ylim = c(0,max(k* yupdf)),
		xlab = (expression(theta)),
		ylab ="Density",
		cex.main=1.5,
		cex.axis=1.5,
		cex.lab=1.5)

	lines(xlpdf, ylpdf,
		col=1,
		lty="solid",
		lwd=2)

	lines(xupdf, k * yupdf,    # times k!
		col=1,
		lty="solid",
		lwd=2)

	if(plot.stat.values){
		meanl   <- MEAN(  drc$dist.lower)
		sdl     <- SD(    drc$dist.lower)
		medianl <- MEDIAN(drc$dist.lower)
		model   <- MODE(  drc$dist.lower)
		abline(v = meanl,         col=2, lty=2, lwd=1)
		abline(v = meanl + sdl/2, col=3, lty=3, lwd=1)
		abline(v = meanl - sdl/2, col=3, lty=3, lwd=1)
		abline(v = model,         col=4, lty=4, lwd=1)
		abline(v = medianl,       col=5, lty=5, lwd=1)
		legend("right", c(paste(c("mean","stdev","mode","median"),
			formatC(c(meanl,sdl,model,medianl)))),
			col = 2:5,
			lty=2:5,
			lwd=2,
			title="values of lower dist",inset = .05)
		meanu   <- MEAN(  drc$dist.upper)
		sdu     <- SD(    drc$dist.upper)
		medianu <- MEDIAN(drc$dist.upper)
		modeu   <- MODE(  drc$dist.upper)
		abline(v = meanu,		  col=6, lty=6, lwd=1)
		abline(v = meanu + sdu/2, col=7, lty=7, lwd=1)
		abline(v = meanu - sdu/2, col=7, lty=7, lwd=1)
		abline(v = modeu, 	      col=8, lty=8, lwd=1)
		abline(v = medianu,       col=9, lty=9, lwd=1)
		legend("topright", c(paste(c("mean","stdev","mode","median"),
			formatC(c(meanu,sdu,modeu,medianu)))),
			col = 6:9,
			lty=6:9,
			lwd=2,
			title="values of upper dist",inset = .05)
	}
	if(makePDF) {dev.off()}
}


Fl.drclass <- function(drc, x, ...)
{
	l <- drc$dist.lower
	u <- drc$dist.upper	
	k <- Kappa(drc)
	Fl<-CDF(l,x) / ( CDF(l,x) + k*(1-CDF(u,x)) )
	return(Fl)
}

FlInv.drclass<-function(drc, p, ...)
{
	search.interval <- c(
						 CDFinv(drc$dist.upper,0.0001),
						 CDFinv(drc$dist.upper,0.9999)
						 )
	FlInv <- uniroot(
	function(x)
		{
			Fl( drc,x ) - p
		},
		lower=min(search.interval),
		upper=max(search.interval),
		f.lower=-1,
		f.upper=1
	)
	return(FlInv$root)
}

Fu.drclass <- function(drc,x, ...)
{
	l <- drc$dist.lower
	u <- drc$dist.upper
	k <- Kappa(drc)
	Fu<-k * CDF(u,x) / ( k * CDF(u,x) + (1 - CDF(l,x)) )
	return(Fu)
}

FuInv.drclass<-function(drc, p, ...)
{
	search.interval <- c(
						 CDFinv(drc$dist.upper,0.0001),
	                     CDFinv(drc$dist.upper,0.9999)
	                     )
	FuInv<-uniroot(
		function(x)
		{
			Fu( drc, x ) - p
		},
		lower=min(search.interval),
		upper=max(search.interval),
		f.lower=-1,
		f.upper=1
	)
	return(FuInv$root)
}

metric.ci.drclass <- function(drc, alpha = 0.05, ...)
{
	i1_alpha <- c(NA,NA)
	i1_alpha[1] <- FuInv(drc, alpha/2)
	i1_alpha[2] <- FlInv(drc, 1-(alpha/2))
	return(i1_alpha)
} # end function metric.ci

metric.width.drclass <- function(drc, alpha = 0.05, ...)
{
	A <- FlInv(drc,(alpha/2)) - FuInv(drc,(alpha/2))
	B <- FlInv(drc,1-(alpha/2)) - FuInv(drc,1-(alpha/2))
	C <- metric.ci(drc,alpha)
	C <- as.numeric(C[2]-C[1])
	relwidth <- (A+B)/(2*C)
	return(relwidth)
} # end function metric.width

metric.shape.drclass <- function(drc, alpha = 0.05, ...)
{
	maxDeltaQ <- -optimize(
						function(x) FuInv(drc,x) - FlInv(drc,x),
						lower=alpha, upper=1-alpha
						)$objective
	C <- metric.ci(drc,alpha)
	C <- as.numeric(C[2]-C[1])
	relshape <- maxDeltaQ/C
	return(relshape)
} # end function rShape1_alpha

metric.mode.drclass <- function(drc, alpha = 0.05, ...)
{
	l  <- drc$dist.lower
	u  <- drc$dist.upper
	k  <- Kappa(drc)
	Mode <- as.numeric(MODE(l))
	ModeHight <- PDF(l,Mode)
	if(!is.na(Mode))
	{
		candidateL <- uniroot(
						function (x) k*PDF(u,x)-ModeHight,
						c(FuInv.drclass(drc,alpha),Mode),
						tol = 0.0001
						)$root
		candidateR <- uniroot(
						function (x) k*PDF(u,x)-ModeHight,
						c(Mode,FlInv.drclass(drc,1-alpha)),
						tol = 0.0001
						)$root
		iMode1_alpha <- candidateR - candidateL
		C <- metric.ci(drc,alpha)
		C <- as.numeric(C[2]-C[1])
		metric.mode <-iMode1_alpha/C
		return(metric.mode)
	} # end !is.na(Mode)
	else return(NA)
} # end function metric.mode

Kappa.drclass<-function(drc, ...)
{
   Kappa <- calc.k(drc$p, drc$q, drc$dist.lower, drc$dist.upper, drc$dist.lower$par, drc$dist.upper$par)$Kappa
   return(Kappa)
}

Lambda.drclass<-function(drc, ...)
{
   Lambda <- calc.k(drc$p, drc$q, drc$dist.lower, drc$dist.upper, drc$dist.lower$par, drc$dist.upper$par)$Lambda
   return(Lambda)
}

metrics.drclass <- function(drc, alpha = 0.05, ...)
{
   mci <-metric.ci(drc, alpha)
   mw <- metric.width(drc, alpha)
   ms <- metric.shape(drc, alpha)
   mm <- metric.mode(drc, alpha)
   return(
   		list(	metric.ci = mci,
   				metric.width=mw,
   				metric.shape = ms,
   				metric.mode =mm
   			)
   	)
}

Fl           <- function(drc,x, ... )  UseMethod("Fl")
Fu           <- function(drc,x, ... )  UseMethod("Fu")
FlInv        <- function(drc,p, ... )  UseMethod("FlInv")
FuInv        <- function(drc,p, ... )  UseMethod("FuInv")
metric.ci    <- function(drc, alpha, ... )  UseMethod("metric.ci")
metric.width <- function(drc, alpha, ... )  UseMethod("metric.width")
metric.shape <- function(drc, alpha, ... )  UseMethod("metric.shape")
metric.mode  <- function(drc, alpha, ... )  UseMethod("metric.mode")
metrics      <- function(drc, alpha, ...)   UseMethod("metrics")
Kappa        <- function(drc,...)  UseMethod("Kappa")
Lambda       <- function(drc,...) UseMethod("Lambda")

################################################################################

calc.k<-function(p = NA, q = NA, dist.lower, dist.upper, par.lower, par.upper)
{
	# overall interval [x.min, x.max]
   	x.min.lower <- CDFinv(dist.lower, 0.001, par.lower)
   	x.max.lower <- CDFinv(dist.lower, 0.999, par.lower)
   	x.min.upper <- CDFinv(dist.upper, 0.001, par.upper)
   	x.max.upper <- CDFinv(dist.upper, 0.999, par.upper) 
	x.min <- min(x.min.lower, x.min.upper)
	x.max <- max(x.max.lower, x.max.upper)
	x <- seq(x.min, x.max, length = 1000000)
	Lambda<-max(PDF(dist.lower, x, par.lower)/PDF(dist.upper, x, par.upper))
	# calculation of k
	k1 <- (p*(1-CDF(dist.lower, q, par.lower)))/(CDF(dist.upper, q, par.upper)*(1-p))
	k2<- (CDF(dist.lower, q, par.lower)*(1-p))/(p*(1-CDF(dist.upper, q, par.upper)))
	Kappa<-max(Lambda, k1, k2)
	if(is.na(p) && is.na(q)){Kappa <- Lambda}
	return(list(Kappa = Kappa, Lambda = Lambda))
} # end function calc.k

###############################################################################

# trans.from.R.to.interval
# =========================

# purpose:
# transforms the real axis to an interval (default: unit interval)

# arguments:
# x:       data to be transformed
# min:     minimum of the interval (default: 0)
# max:     maximum of the interval (default: 1)

# output:
# transformed data

trans.from.R.to.interval <- function(x, min = 0, max = 1)
{
	y <- 0.5*(min+max) + (max-min)/pi*atan(x)
	return(y)
}

################################################################################

# trans.from.interval.to.R
# ===========================

# purpose:
# transforms an interval (default: unit interval) to the real axis

# arguments:
# y:       data to be transformed
# min:     minimum of the interval (default: 0)
# max:     maximum of the interval (default: 1)

# output:
# transformed data

trans.from.interval.to.R <- function(y, min = 0, max = 1)
{
	x <- tan(0.5*pi*(2*y-max-min)/(max-min))
	return(x)
}

################################################################################

# trans.from.R.to.Rplus
# ===========================

# purpose:
# transforms the real axis (default: unit interval) to the positive real axis

# arguments:
# x:       data to be transformed
# min:     [min,Inf] (default: [0,Inf])

# output:
# transformed data

trans.from.R.to.Rplus <- function(x, min = 0)
{
	y <- exp(x)+min
	return(y)
}

################################################################################

# trans.from.Rplus.to.R
# ===========================

# purpose:
# transforms the positive real axis (default: [0,Inf]) to the real axis

# arguments:
# y:       data to be transformed
# min:     [min,Inf] (default: [0,Inf])

# output:
# transformed data

trans.from.Rplus.to.R <- function(y, min = 0)
{
	x <- log(y-min)
	return(x)
}

################################################################################

# trans.from.interval.to.interval
# =================================

# purpose:
# transforms an interval ([min1,max1]) to an interval ([min2,max2]) 

# arguments:
# x:       data to be transformed
# min:     minimum of the interval (default: 0)
# max:     maximum of the interval (default: 1)

# output:
# transformed data

trans.from.interval.to.interval <- function(x, min1 = 0, max1 = 1, min2 = 0, max2 = 1)
{
	y <- (x-min1) * (max2-min2)/(max1-min1) + min2
	return(y)
}

################################################################################
aberr.l.bfgs.b <- function( par.optim, p, q, dist.lower, dist.upper, start.dist.lower.par, start.dist.upper.par, ...)
{	
    if( !is.na( start.dist.lower.par )[1] )
	{	
		par.lower <- start.dist.lower.par # ini
		for(i in 1:length(start.dist.lower.par))
		{	
			eval(parse(text = paste("par.lower[",i,"]<-par.optim[",i,"]",sep = "")))
		}
		z <- length(start.dist.lower.par)
	}
	if( is.na( start.dist.lower.par )[1] )
	{
		par.lower <- dist.lower$par
		z <- 0
	}
	if( !is.na( start.dist.upper.par )[1]  )
	{
		par.upper <- start.dist.upper.par # ini
		for(i in 1:length( start.dist.upper.par ))
		{	
			eval(parse(text = paste("par.upper[",i,"]<-par.optim[",z+i,"]",sep="")))
		}
	}
	if( is.na( start.dist.upper.par )[1]  )
	{
		par.upper <- dist.upper$par
	}
	par.lower <- mergePar( par.lower, dist.lower$par, dist.lower$par.names )
	par.upper <- mergePar( par.upper, dist.upper$par, dist.upper$par.names )

	k<-calc.k( p, q, dist.lower, dist.upper, par.lower, par.upper )$Kappa
	cat(paste("\n Kappa:",k,"\n\n"))
	cat(paste(par.lower,"\n"))
	cat(paste(par.upper,"\n"))
	if( k < 999 ){return(as.numeric(k))}
	if( k >= 999 ) 
	{
		return( 999 )
	}
	
} # end function aberr.l.bfgs.b

################################################################################
aberr.nelder.mead <- function( par.optim, p, q, dist.lower, dist.upper, start.dist.lower.par, start.dist.upper.par, M, ...)
{	
	for( i in 1:length(par.optim))
	{
		if(M[1,i] == -Inf)
		{
			if(!M[2,i] == Inf)
			{
				par.optim[i] <- -trans.from.R.to.Rplus(par.optim[i], -M[2,i])
			}
		}
		else
		{
			if(M[2,i] == Inf)
			{
				par.optim[i] <- trans.from.R.to.Rplus(par.optim[i], M[1,i])
			}
			else
			{
				par.optim[i] <- trans.from.R.to.interval(par.optim[i], M[1,i], M[2,i])
			}
		}
	}
	if( !is.na( start.dist.lower.par )[1] )
	{	
		par.lower <- start.dist.lower.par # ini
		for(i in 1:length(start.dist.lower.par))
		{	
			eval(parse(text = paste("par.lower[",i,"]<-par.optim[",i,"]",sep = "")))
		}
		z <- length(start.dist.lower.par)
	}
	if( is.na( start.dist.lower.par )[1] )
	{
		par.lower <- dist.lower$par
		z <- 0
	}
	if( !is.na( start.dist.upper.par )[1]  )
	{
		par.upper <- start.dist.upper.par # ini
		for(i in 1:length( start.dist.upper.par ))
		{	
			eval(parse(text = paste("par.upper[",i,"]<-par.optim[",z+i,"]",sep="")))
		}
	}
	if( is.na( start.dist.upper.par )[1]  )
	{
		par.upper <- dist.upper$par
	}
	par.lower <- mergePar( par.lower, dist.lower$par, dist.lower$par.names )
	par.upper <- mergePar( par.upper, dist.upper$par, dist.upper$par.names )

	k<-calc.k( p, q, dist.lower, dist.upper, par.lower, par.upper )$Kappa
	cat(paste("\n Kappa:", k, "\n\n"))
	cat(paste(par.lower, "\n"))
	cat(paste(par.upper, "\n"))
	if( k < 999 ){return(as.numeric(k))}
	if( k >= 999 ) 
	{
		return( 999 )
	}
	
} # end function aberr.nelder.mead 

################################################################################

process.elidat <- function(p = p, q = q, dist.lower, dist.upper, start.dist.lower.par = NA, start.dist.upper.par = NA, ...)
{
	if( is.na(start.dist.lower.par) && is.na(start.dist.upper.par) )
	{
		calculated <- calc.k( p, q, dist.lower, dist.upper, dist.lower$par, dist.upper$par )
	    drc                      <- list()
	    drc$name                 <- "Non-Optimised"
	    drc$p                    <- p
	    drc$q                    <- q
	    drc$dist.lower           <- dist.lower
	    drc$dist.upper           <- dist.upper
	    class(drc)               <- "drclass"
		return(
			list(
			    drc = drc,
				k = calculated$Kappa,
				optimised.par.lower = NA,
				optimised.par.upper = NA,
				par.lower = dist.lower$par,
				par.upper = dist.upper$par
			)
		)
	} # end if( is.na(start.dist.lower.par) && is.na(start.dist.upper.par) )
	else
	{
	   par.optim <- as.vector(na.omit(c(start.dist.lower.par,start.dist.upper.par)))
	   id <- as.numeric(
	             pmatch( names(c(start.dist.lower.par, start.dist.upper.par )),
	                     names(c(dist.lower$par, dist.upper$par )))
	          )
	   lower <- rep(NA,length(id))
	   upper <- rep(NA,length(id))
	   z <- 1       
	   for(i in c(id)){ 
          lower[z]     <- c(dist.lower$par.range[,1],dist.upper$par.range[,1])[i]
		  upper[z]     <- c(dist.lower$par.range[,2],dist.upper$par.range[,2])[i]
		  z <- z+1
		} # end for i
	    M <- matrix(c(lower,
                      upper),
                      byrow = TRUE,
                      ncol = length(par.optim)
                    )
	    for( z in 1:10 )
	    {
		    res<-optim(	par.optim,
					aberr.l.bfgs.b,
					gr = NULL,
					p, q, dist.lower, dist.upper, start.dist.lower.par, start.dist.upper.par, M,
					method = "L-BFGS-B",
					lower = lower,
					upper = upper,
					control=c(factr=1e7,maxit=10000)
					)
		    for(i in 1:length(par.optim))
		    {
		    	par.optim[i] <-res$par[i] #*(1+rnorm(1,0,0.000001))
		    }# end i
		    # whatever to R
		    for( i in 1:length(par.optim))
		    {
			    if(M[1,i] == -Inf)
			    {
				    if(!M[2,i] == Inf)
				    {
					    par.optim[i] <- trans.from.Rplus.to.R(-par.optim[i],-M[2,i])
			    	}
			    }
			    else
			    {
				    if(M[2,i] == Inf)
				    {
					    par.optim[i] <- trans.from.Rplus.to.R(par.optim[i],M[1,i])
			    	}
				    else
				    {
					    par.optim[i] <- trans.from.interval.to.R(par.optim[i],M[1,i],M[2,i])
				    }
			    }
		    }
		    res<-optim(	par.optim,
					aberr.nelder.mead,
					gr = NULL,
					p, q, dist.lower, dist.upper, start.dist.lower.par, start.dist.upper.par, M,
					method = "Nelder-Mead",
					control=c(maxit=2000,abstol=0.000001)
					)
		    for(i in 1:length(par.optim))
		    {
			    par.optim[i] <-res$par[i] #*(1+rnorm(1,0,0.000001))
	  	    }# end i
	  	    for( i in 1:length(par.optim))
		    {
			    if(M[1,i] == -Inf)
			    {
				    if(!M[2,i] == Inf)
				    {
					    par.optim[i] <- -trans.from.R.to.Rplus(par.optim[i],-M[2,i])
				    }
			    }
			    else
			    {
				    if(M[2,i] == Inf)
				    {
					par.optim[i] <- trans.from.R.to.Rplus(par.optim[i],M[1,i])
				    }
				    else
				    {
					    par.optim[i] <- trans.from.R.to.interval(par.optim[i],M[1,i],M[2,i])
				    }
			    }
		    } # end i	
		    print(paste("The optimisation algorithm was ",z," time(s) repeated."))
		    if(res$convergence==0)
		    {
			    break
	    	}
	    } # end z
	    if( !is.na( start.dist.lower.par )[1] )
	    {
		    optimised.par.lower <- start.dist.lower.par # ini
		    for(i in 1:length( start.dist.lower.par ))
		    {	
			    eval(parse(text = paste("optimised.par.lower[",i,"]<-par.optim[",i,"]",sep = "")))
    		}
		    par.lower <- mergePar(optimised.par.lower, dist.lower$par, dist.lower$par.names)
		    z <- length(start.dist.lower.par)
	    }
	    if( is.na( start.dist.lower.par )[1] )
	    {
		    z <- 0
		    optimised.par.lower <- NA
		    par.lower           <- dist.lower$par
		    #names(par.lower)    <- dist.lower$par.names
	    }	
	    if( !is.na( start.dist.upper.par )[1] )
	    {
		    optimised.par.upper <- start.dist.upper.par # ini
		    for(i in 1:length( start.dist.upper.par ))
		    {
			    eval(parse(text = paste("optimised.par.upper[",i,"]<-par.optim[",z+i,"]",sep="")))
		    }
	        par.upper <- mergePar(optimised.par.upper, dist.upper$par, dist.upper$par.names)
	    }
	    if( is.na( start.dist.upper.par )[1] )
	    {
		    optimised.par.upper <- NA
		    par.upper           <- dist.upper$par
		    #names(par.upper)    <- dist.upper$par.names
	    }
	    dist.lower$par       <- par.lower
	    dist.upper$par       <- par.upper
	    # fill DRC
	    drc                      <- list()
	    drc$name                 <- "Optimised"
	    drc$range                <- dist.upper$range(dist.upper$par)
	    drc$p                    <- p
	    drc$q                    <- q
	    drc$dist.lower           <- dist.lower
	    drc$dist.upper           <- dist.upper
	    class(drc)               <- "drclass"
	    return(
		    list
		    (	
			    drc = drc,
		    	optimised.par.lower  = optimised.par.lower,
			    optimised.par.upper  = optimised.par.upper
		    )
	    )
	} # end else

} # end function process.elidat

################################################################################

# Tests
#trans.arctan <- trans.arctan.create(c(0,10))
#dist.uniform <- dist.uniform.create(c(0,10))
#dist.lower  <- dist.trans.create(dist.lower,trans.tan)
#dist.weibull  <- dist.weibull.create(c(1,200))
#dist.student <- dist.student.create(c(0,1,1000))
#dist.upper  <- dist.trans.create(dist.student,trans.tan)
#trans.tan <- trans.tan.create(c(33,99))
#dist.logistic <- dist.logistic.create(c(0, 0.1))
#dist.lower  <- dist.trans.create(dist.logistic,trans.arctan)
#x11()
#plot(dist.lower,what="CDFinv",plot=FALSE)
#x11()
#plot(dist.lower,par=c(0,1,NA))
#x11()
#plot(dist.upper,par=c(0,1,1000,0,150))
#x11()
#plot(dist.upper)

