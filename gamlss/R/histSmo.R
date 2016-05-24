# needs this to read the data reading the data 
# Paul Eilers' density estimation
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
histSmo <- function(y, 
                lambda = NULL, 
                    df = NULL, 
                 order = 3, 
                 lower = NULL,  
                 upper = NULL, 
                  type = c("freq", "prob"),  
                  plot = FALSE, 
                breaks = NULL,  
              discrete = FALSE, ...)
{
histSmocall <- match.call() 
       type <- match.arg(type) 
     if (is.null(lambda)&&is.null(df))# estimate lambda
     {
m1 <- histSmoP(y=y, lambda=lambda, df=df, order=order, lower=lower, upper=upper, type=type, plot=plot, breaks = breaks, discrete = discrete, ...) 
     }
       if (!is.null(lambda))# specify lambda
     {
m1 <- histSmoO(y=y, lambda=lambda,order=order, lower=lower, upper=upper, type=type, plot=plot, breaks = breaks, discrete = discrete, ...)  
     }   
       if (!is.null(df))# specify df's
     {
m1 <- histSmoC(y=y, df=df, lower=lower, upper=upper, type=type, plot=plot, breaks = breaks, discrete = discrete, ...)  
     } 
 m1$call <- histSmocall
invisible(m1)
}
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# the original Paul's function
# it works with fixed lambda
# The Poisson smoother
histSmoO <- function(y, 
                    lambda = 1, order=3, 
                     lower = NULL, 
                     upper = NULL, 
                      type = c("freq", "prob"), 
                      plot = FALSE, 
                    breaks = NULL,
                  discrete = FALSE,
                     ...)
{
# local functions	
  histsm <- function(y, lambda = 10, d = 2)
 {
  # Penalty stuff
  m <- length(y)
  E <- diag(m)
  D <- diff(E, diff = d)            
  G <- lambda * t(D) %*% D  
  # Initialie
  eta <- log(y + 0.5)
  mu0 <- 0
  # Iterate
  for (it in 1:20) 
    {
     mu <- exp(eta)
    dmu <- max(abs(mu - mu0))
    mu0 <- mu
     W <- diag(mu)
    eta  <- solve(W + G, y - mu + mu * eta)
    if (dmu < 1e-5) break
    }
  return(mu)
 }
#--------------------------------
# this function is for descrete variable smooting
# it take a discrete variable and creates values and frequencies 	
getTabs <- function(y, freq=NULL, lower=min(y), upper=max(y))
{
        x1 <- seq(lower, upper, by=1)
	      fy <- if (is.null(freq)) factor(y,levels=x1) else factor(rep(y,freq),levels=x1)
	  tabley <- xtabs(~fy)
	#notresy <- if (is.null(freq)) factor(y) else factor(rep(y,freq)) 
       out <- data.frame(x=x1, freq=as.numeric(tabley))# get the table    
out
}  
#--------------------------------
# main function starts here 	
   histSmocall <- match.call() 
          type <- match.arg(type)
            Ry <- range(y) #ok
            ly <- length(y) # ok
            Ey <- (Ry[2]-Ry[1])*.10 #ok
     lower.lim <- if (is.null(lower)) (Ry[1]-Ey) else lower #ok
     upper.lim <- if (is.null(upper))  (Ry[2]+Ey) else upper #ok
        breaks <- if (is.null(breaks)) floor(length(y)/10) else breaks#ok
  if (discrete)
        {
        	lower <- if (is.null(lower)) min(y) else lower
        	upper <- if (is.null(upper)) max(y)+1 else upper
              hst <- getTabs(y, lower=lower, upper=upper)
                x <- as.numeric(hst$x)
                Y <- hst$freq    
        }
        else
        {
              hst <- hist(y, breaks = seq(lower.lim,upper.lim , length = breaks), plot=FALSE)
                x <- hst$mids
               dx <- (x[2] - x[1])
                Y <- hst$counts
        }  
          # hst <- hist(y, breaks = seq(lower.lim,upper.lim , length = breaks), plot=FALSE)
           #  x <- hst$mids
            #dx <- (x[2] - x[1])
            # Y <- hst$counts
            mu <- histsm(Y, lambda, order)
       ncounts <- sum(Y)
         cdf <- cumsum(mu) / sum(mu)
      cdfFun <- stepfun(x, c(0,cdf))
        invcdfFun <- if (discrete) stepfun(cdf, c(x, max(x)))  else splinefun(cdf,x)
if (plot)
 {
  if (discrete)
   {
  switch(type, "freq"={ 
  	     r<-barplot(Y, fg = "blue", col="gray", axis.lty=1, 
            border="blue", col.axis = "blue4",col.main = "blue4", 
            col.lab = "blue4",  ylab="Frequency", names.arg=as.character(x),
            ...)
            lines(r, mu, col = 'red', lty = 1, lwd = 1, type="h")
            points(r, mu, col="red")
                      },
                  "prob"= { 
        r <- barplot(Y/sum(Y), fg = "blue", col="gray", axis.lty=1, 
             border="blue", col.axis = "blue4",col.main = "blue4", 
            col.lab = "blue4",   ylab="Density", names.arg=as.character(x),
                ...)
            lines(r, (mu/sum(mu)), col = 'red', lty = 1, lwd = 3)
            points(r, (mu/sum(mu)), col="red")
                      }) 		
   }	
  else
   {
  switch(type, "freq"={plot(hst, ...);lines(x, mu, col = 'red', lty = 1, lwd = 3)},
                "prob"= {plot(hst, freq=FALSE,...); lines(x, (mu/sum(mu))/dx, col = 'red', lty = 1, lwd = 3)}) 	
   }
 }

         density <- if (discrete) (mu/sum(mu)) else (mu/sum(mu))/dx	
             out <- list(x=x, counts=Y, density=density, hist=hst, cdf=cdfFun, invcdf=invcdfFun, call= histSmocall, discrete=discrete)
      class(out) <- "histSmo"
  invisible(out) 

}    
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# The Poisson smoother using cubic splines and GAMLSS
# it works but it needs df's'
#histSmoC <- function(y, df=10, lower=NULL, upper=NULL, type=c("freq", "prob"),  save=FALSE, plot=TRUE, breaks=NULL)
histSmoC <- function(y,  
                      df = 10,  
                   lower = NULL, 
                   upper = NULL, 
                    type = c("freq", "prob"), 
                    plot = FALSE, 
                  breaks = NULL, 
                discrete = FALSE,
                     ...)

{
#-----------------------------
# local function	
# this function is for descrete variable smooting
# it take a descrete variable and creates values and frequencies 	
getTabs <- function(y, freq=NULL, lower=min(y), upper=max(y))
{
          x1 <- seq(lower, upper, by=1)
	      fy <- if (is.null(freq)) factor(y,levels=x1) else factor(rep(y,freq),levels=x1)
	  tabley <- xtabs(~fy)
	#notresy <- if (is.null(freq)) factor(y) else factor(rep(y,freq)) 
         out <- data.frame(x=x1, freq=as.numeric(tabley))# get the table    
out
}  
#------------------------------
# main function starts here 	
 histSmocall <- match.call() 
#        require(gamlss)    
          type <- match.arg(type)
            Ry <- range(y)
            ly <- length(y) 
            Ey <- (Ry[2]-Ry[1])*.10
     lower.lim <- if (is.null(lower)) (Ry[1]-Ey) else lower
     upper.lim <- if (is.null(upper))  (Ry[2]+Ey) else upper    
        breaks <- if (is.null(breaks)) max(floor(length(y)/10), 101) 
                 else breaks
      # 
        if (discrete)
        {
        	lower <- if (is.null(lower)) min(y) else lower
        	upper <- if (is.null(upper)) max(y)+1 else upper
              hst <- getTabs(y, lower=lower, upper=upper)
                x <- as.numeric(hst$x)
                Y <- hst$freq    
        }
        else
        {
              hst <- hist(y, breaks = seq(lower.lim,upper.lim , length = breaks), plot=FALSE)
                x <- hst$mids
               dx <- (x[2] - x[1])
                Y <- hst$counts
        }  
          ncounts <- sum(Y)
               m1 <- gamlss(Y~cs(x, df=df), trace=FALSE, family=PO)
               mu <- fitted(m1)
              cdf <- cumsum(mu) / sum(mu)
           cdfFun <- stepfun(x, c(0,cdf))
        invcdfFun <- if (discrete) stepfun(cdf, c(x, max(x)))  else splinefun(cdf,x)
if (plot)
 {
  if (discrete)
   {
  switch(type, 
         "freq"={ r<-barplot(Y, fg = "blue", col="gray", axis.lty=1, 
          border="blue", col.axis = "blue4",col.main = "blue4", 
          col.lab = "blue4",  ylab="Frequency", names.arg=as.character(x), ...)
          lines(r, mu, col = 'red', lty = 1, lwd = 1, type="h")
                                    points(r, mu, col="red")
                },
          "prob"= { r <- barplot(Y/sum(Y), fg = "blue", col="gray", axis.lty=1, 
              border="blue", col.axis = "blue4",col.main = "blue4", 
            col.lab = "blue4",   ylab="Density", names.arg=as.character(x), ...)
            lines(r, (mu/sum(mu)), col = 'red', lty = 1, lwd = 3)
                                   points(r, (mu/sum(mu)), col="red")
                                    }) 		
   }	
  else
   {
  switch(type, "freq"={plot(hst, ...);lines(x, mu, col = 'red', lty = 1, lwd = 3)},
                "prob"= {plot(hst, freq=FALSE, ...); lines(x, (mu/sum(mu))/dx, col = 'red', lty = 1, lwd = 3)}) 	
   }
 }
  {
         density <- if (discrete) (mu/sum(mu)) else (mu/sum(mu))/dx	
             out <- list(x=x, counts=Y, density=density, hist=hst, cdf=cdfFun, invcdf=invcdfFun, model=m1, call= histSmocall, discrete=discrete)
      class(out) <- "histSmo"
      invisible(out) 
  }       
}
#---------------------------------------------------------------#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# The Poisson smoother using P-splines and GAMLSS
# it works but breaks need to be more that 100 otherwise the knots=10 of the pb() ia take off
# one solution is to change the pb( function condition from length(Y)<100 to <99)
# this is a new version of histSmpoP() allowing descete variable    
# 3-Dec-2011 Mikis  
#-------------------------------------------------------------  
histSmoP <- function(y, 
                  lambda = NULL, 
                      df = NULL,  
                   order = 3, 
                   lower = NULL, 
                   upper = NULL, 
                    type = c("freq", "prob"), 
                    plot = FALSE, 
                  breaks = NULL, 
                discrete = FALSE, ...)
{
#-----------------------------
# local function	
# this function is for descrete variable smooting
# it take a descrete variable and creates values and frequencies 	
getTabs <- function(y, freq=NULL, lower=min(y), upper=max(y))
{
        x1 <- seq(lower, upper, by=1)
	      fy <- if (is.null(freq)) factor(y,levels=x1) else factor(rep(y,freq),levels=x1)
	  tabley <- xtabs(~fy)
	#notresy <- if (is.null(freq)) factor(y) else factor(rep(y,freq)) 
       out <- data.frame(x=x1, freq=as.numeric(tabley))# get the table    
out
}  
#------------------------------
# main function starts here 	
#        require(gamlss)    
   histSmocall <- match.call() 
          type <- match.arg(type)
            Ry <- range(y)
            ly <- length(y)
            if (ly<100) warning("with less that 100 observations reduce the default argument breaks which is 101")  
            Ey <- (Ry[2]-Ry[1])*.10
     lower.lim <- if (is.null(lower)) (Ry[1]-Ey) else lower
     upper.lim <- if (is.null(upper))  (Ry[2]+Ey) else upper    
        breaks <- if (is.null(breaks)) max(floor(length(y)/10), 101) 
                 else breaks
      # 
        if (discrete)
        {
        	lower <- if (is.null(lower)) min(y) else lower
        	upper <- if (is.null(upper)) max(y)+1 else upper
            hst <- getTabs(y, lower=lower, upper=upper)
              x <- as.numeric(hst$x)
              Y <- hst$freq    
        }
        else
        {
              hst <- hist(y, breaks = seq(lower.lim,upper.lim , length = breaks), plot=FALSE)
                x <- hst$mids
               dx <- (x[2] - x[1])
                Y <- hst$counts
        }  
          ncounts <- sum(Y)
               m1 <- gamlss(Y~pb(x, df=df, lambda=lambda, inter=50), trace=FALSE, family=PO)
               mu <- fitted(m1)
              cdf <- cumsum(mu) / sum(mu)
           cdfFun <- stepfun(x, c(0,cdf))
        invcdfFun <- if (discrete) stepfun(cdf, c(x, max(x)))  else splinefun(cdf,x)
if (plot)
 {
  if (discrete)
   {
  switch(type, "freq"={ 
  	     r<-barplot(Y, fg = "blue", col="gray", axis.lty=1, 
            border="blue", col.axis = "blue4",col.main = "blue4", 
            col.lab = "blue4",  ylab="Frequency", names.arg=as.character(x),
            ...)
          lines(r, mu, col = 'red', lty = 1, lwd = 1, type="h")
          points(r, mu, col="red")
                        },
                "prob"= { 
        r <- barplot(Y/sum(Y), fg = "blue", col="gray", axis.lty=1, 
             border="blue", col.axis = "blue4",col.main = "blue4", 
            col.lab = "blue4",   ylab="Density", names.arg=as.character(x),
                                    ...)
            lines(r, (mu/sum(mu)), col = 'red', lty = 1, lwd = 3)
            points(r, (mu/sum(mu)), col="red")
                                    }) 		
   }	
  else
   {
  switch(type, "freq"={plot(hst, ...);lines(x, mu, col = 'red', lty = 1, lwd = 3)},
                "prob"= {plot(hst, freq=FALSE, ...); lines(x, (mu/sum(mu))/dx, col = 'red', lty = 1, lwd = 3)}) 	
   }
 }
 
         density <- if (discrete) (mu/sum(mu)) else (mu/sum(mu))/dx	
             out <- list(x=x, counts=Y, density=density, hist=hst, cdf=cdfFun, invcdf=invcdfFun, model=m1, call= histSmocall, discrete=discrete)
      class(out) <- "histSmo"
      invisible(out) 
         
}
#---------------------------------------------------------------
#---------------------------------------------------------------
plot.histSmo <- function(x, type=c("hist", "cdf", "invcdf"), ...) 
{
	     type <- match.arg(type)
	 discrete <- x$discrete
	if (discrete)
   {
   		switch(type,"hist"={r<-barplot(x$counts/(sum(x$counts)), fg = "blue", col="gray", axis.lty=1,  border="blue", col.axis = "blue4",col.main = "blue4", 
  col.lab = "blue4",  ylab="Frequency", names.arg=as.character(x$x))
          lines(r, x$density, col = 'red', lty = 1, lwd = 1, type="h")
          points(r, x$density, col="red")}, 
	                "cdf"={plot(x$cdf, ylab="cdf")},
	            "invcdf"={plot(x$invcdf, ylab="invcdf", xlab="p")}
	        ) 
   }
   else
   {
   	switch(type,"hist"={plot(x$hist, freq=FALSE); lines(x$x, x$density, col = 'red', lty = 1, lwd = 3)}, 
	               "cdf"={plot(x$cdf, ylab="cdf")},
	            "invcdf"={plot(x$invcdf, 0, 1, ylab="invcdf", xlab="p")}
	        )   
    } 	 
}
#-----------------------------------------------------------------------------
#------------------------------------------------------------------------------
lines.histSmo <- function(x, type=c("hist", "cdf", "invcdf"), ...) 
{
  type <- match.arg(type)
  discrete <- x$discrete
  if (discrete)
  {
    stop("lines is only for continuous data fitted histSmo objects")
#     r<-barplot(x$counts/(sum(x$counts)))
#     switch(type,"hist"={lines(r, x$density, col = 'red',  ...)},
#                 "cdf"={lines(x$cdf, ylab="cdf")},
#              "invcdf"={lines(x$invcdf, ylab="invcdf", xlab="p")}) 
  }
  else
  {
    switch(type,"hist"={lines(x$x, x$density, ...)}, 
           "cdf"={lines(x$cdf, ylab="cdf")},
           "invcdf"={lines(x$invcdf, 0, 1, ylab="invcdf", xlab="p")}
    )   
  } 	 
}

#-------------------------------------------------------------------------------
print.histSmo <- function (x, digits = NULL, ...) 
{
  cat("\nCall:\n\t", deparse(x$call), "\n\n", sep = "")
      #"\n\nData: ", x$data.name, 
 #    " (", length, " obs.);", "\tBandwidth 'bw' = ", formatC(x$bw, 
#                                                           digits = digits), 
 # print(summary(as.data.frame(x[c("x", "y")])), digits = digits, 
 #       ...)
  invisible(x)
}
#-------------------------------------------------------------------------------
#Data = read.csv(path1,header = T)   
#x11()

#for (j in 1:nrow(Data))  out[j] <-np.pdf(j)
