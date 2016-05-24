vitalsim <-function(vrmeans, vrvars, corrin, corrout, makemx, n0, yrspan, Ne=500, tmax=50, runs=500, vrtypes=NULL, vrmins=NULL, vrmaxs=NULL,sumweight=NULL)
{
x1<-length(vrmeans)  
x2<-length(n0)
x3<-dim(corrin)[1]
## calcualte np and np2 
np<-x3
np2<-x1-x3
np3 <- np*yrspan   ## CORRECTION in email from Doak 8/4/07

## create some reasonable defaults if not specified
 if(missing(vrmins)){vrmins=rep(0,x1)}
 if(missing(vrmaxs)){vrmaxs=rep(0,x1)}
 if(missing(vrtypes)){vrtypes=rep(1,x1)} # 1 = beta, 2 = stretched beta, 3 = lognormal
 if(missing(sumweight)){sumweight=rep(1,x2)}
## some error checking -- could add more
  if(length(vrtypes)!=x1){ stop("vrtypes is not the same length as vrmeans!", call.=FALSE)}
  if(length(vrmins)!=x1){ stop("vrmins is not the same length as vrmeans!", call.=FALSE)}
  if(length(vrmaxs)!=x1){ stop("vrmaxs is not the same length as vrmeans!", call.=FALSE)}
  if(length(sumweight)!=x2){ stop("sumweight is not the same length as n0!", call.=FALSE)}
vrs <- vrmeans	 
meanmx<-makemx(vrs)                      # create mean matrix using function
#
Nstart <- sum(n0)	                 # starting population number
lam0 <- lambda(meanmx)                   # find the deterministic population growth rate (mean matrix)
#-------------------------------------------------------------------------------
# this section makes sets of beta or strecthed beta values to choose
# from during the simulations. It makes 99 values for 1
# increments of Fx for each parameter 
  print('Generating beta distributed values for vital rates', quote=FALSE)
  parabetas <- matrix(0,99,np+np2)
  for (ii in 1:(np+np2))
  {
     if (vrtypes[ii] != 3)
     {
   	 for (fx99 in 1:99)
         {
    	     if (vrtypes[ii] == 1)
             {
  		parabetas[fx99,ii] <- betaval(vrmeans[ii],sqrt(vrvars[ii]),fx99/100)
             }
  	     if (vrtypes[ii] == 2)
             {
  		 parabetas[fx99,ii] <- stretchbetaval(vrmeans[ii],sqrt(vrvars[ii]),
                                                      vrmins[ii], vrmaxs[ii], fx99/100)
             }
   	 } 
     } 
  } 
#--------creating and using the big correlation matrix, M--------
# this set of loops makes the big correlation matrix (M)
# with multi-year correlations: the if statements are used to
# estimate the correct correlations with increasing time lags,
# always assuming that all long-time-lag correlations are only
# caused by within-year and one-time-step correlations
#
# need function for simple matrix powers 
"%^%" <- function(mat, pow)
{
   result <- diag(nrow(mat))
   while (pow > 0)
   {
     result <- result %*% mat
     pow <- pow - 1
   }
   result
}

  print('Calculating the multi-year correlation matrix', quote=FALSE)


  M <- matrix(,yrspan*np,yrspan*np)  # initialize the big correlation matrix (M)
  for (ii in 1:yrspan)
  {
     for (jj in 1:yrspan)
     {
    	if (ii == jj)
        {
           litmx <- corrin
        }
        else
        {
           litmx <- corrout
        }
        expo<-1
  	if (ii > jj)
        {
           expo <- ii-jj
           litmx <- litmx %^% expo
        }
  	if (ii < jj)
        {
           expo <- jj-ii
           litmx <- (t(litmx)) %^%expo
        }
	for (ip in 1:np)
        {
   	   for (jp in 1:np)
           {
  	       M[(ip+np*(ii-1)),(jp+np*(jj-1))] <- litmx[ip,jp]
  	   } 
  	} 
     } 
  } 
# get the eigenvalues for calculating the M12 matrix
  ev <- eigen(M)
  d  <- ev$val;  W <- ev$vec           # eigenvalues (diagonal matrix) and vectors
  o <- rev(order(Mod(d)))              # re-order by largest real part to match Matlab??
  d <- d[o]
  W <- W[,o]
## NOTE: this section will not be evaluated using matlab code in Box 8.10
##  since min(d) will always be > 0 if computed using complex input..
  checkeig <- min(d)	              # check for negative eigenvalues

 if (checkeig < 0)
  {
        print('Correcting negative eigenvalues', quote=FALSE)

      maxneg <- abs(min(d[d<0]))      # the largest negative eigenvalue
      d[d <= maxneg] <- 0             # sets negatives and small positive values = 0
      d <- diag(d)
      newfullmx <- W %*% d %*% t(W)   # make a corrected matrix
      for (ii in 1:np3)                # CORRECTION in email 8/4/07
      {                                 ##change from covariances to correlations
  	  for (jj in 1:np3)
          {
  	      if (newfullmx[ii,ii] == 0 | newfullmx[jj,jj] == 0)
              {
                  newfullmx[ii,jj] <- 0
              }
  	      else
              {
  	          newfullmx[ii,jj] <- newfullmx[ii,jj]/((newfullmx[ii,ii] * newfullmx[jj,jj])^0.5)
              }
          } 
      } 
      ev <- eigen(newfullmx)
      d <- ev$val; W <- ev$vec
      o<-rev(order(Mod(d)))
      d<-d[o]
      W<-W[,o]
  }                       

  d <- diag(d)
  M12 <- W %*% (abs(d)^0.5) %*% t(W)  # the M^(1/2) matrix
  sz <- nrow(M12)                     # the total number of lines in M12
  # get the lines from the middle of M12 to use to generate correlations
  startcase <- round(yrspan/2) * np + 1   #
  zvalold  <-  Re(M12[startcase:(startcase + np - 1),])
  zvalnew<-    Re(M12[ (startcase+np):(startcase + 2*np - 1),])
  newns <- matrix(rnorm(sz), ncol=1) #
  oldxy <- zvalold %*% newns #
  #-----end of: creating and using the big correlation matrix------
#
# Runs to get growth rate and extinction risk
print('Running projections to get growth rate and extinction risk', quote=FALSE)
normresults <- c()
PrExt     <- matrix(0,tmax,1)    # the extinction time tracker
logLam    <- matrix(0,runs,1)    # the tracker of log-lambda values
stochLam  <- matrix(0,runs,1)    # tracker of stochastic lambda values
 for (xx in 1:runs)
 {   
    if(xx==1 || xx %% 10 == 0){print(paste("  Starting run", xx), quote=FALSE)}
    nt <- n0 # start at initial population vector
    extinct <- 0
    for (tt in 1:tmax)
    {
     	newns <- matrix(newns[(np+1):sz,],(sz-np),1)  # make random normals
    	newns <- rbind(newns, matrix(rnorm(np),np,1))
     	newxy <- zvalnew %*% newns                    # make new set of correlated normals
        # these lines save normals to check correlations :
        normresults <- rbind(normresults, cbind(t(oldxy), t(newxy)))
        oldxy <- newxy
        # adds in randoms for uncorrelated vital rates.
        yrxy  <- rbind(newxy,matrix(rnorm(np2),np2,1))
        # find vital rate values
        vrs   <- matrix(0,np+np2)            # initialize vrs
        yrxy1 <- matrix(yrxy[vrtypes != 3])  # if not a lognormal rate
        yrxy2 <- matrix(yrxy[vrtypes == 3])  # if lognormal rate
        index <- c()                         # initilize index  
        index <- round(100 * pnorm(yrxy1,0,1,TRUE,FALSE))
        index[index == 0]   <- 1            # round at extremes
        index[index == 100] <- 99
        vrs[vrtypes != 3]   <- diag(parabetas[index,vrtypes != 3])  # find stored value
                                            # calculate a lognormal value
        #vrs[vrtypes == 3]   <- lnorms(vrmeans[vrtypes == 3],vrvars[vrtypes == 3], yrxy2)
         vrs[vrtypes == 3]   <- lnorms(yrxy2, vrmeans[vrtypes == 3],vrvars[vrtypes == 3]) # updated lnorms
        
        mx  <- makemx(vrs)                  # use matrix definition function to make yearly matrix
        nt <- mx %*% nt	 	            # multiply matrix by the population vector
        if (extinct == 0)                   # check for extinction
        {        
            Ntot = sumweight %*% nt	    # compute weighted sum of current densities
            if (Ntot <= Ne)
            {
  	       PrExt[tt] = PrExt[tt] + 1
  	       extinct <- 1
    	    } 
        } 
    } 
    logLam[xx]    <- (1/tmax) * log(sum(nt)/Nstart)   # calculate loglambda
    stochLam[xx]  <- (sum(nt) / Nstart)^(1/tmax)      # and stoch. lambda  Not used??
 } # END runs loop

 CDFExt  <- cumsum(PrExt/runs)               # make the extinction CDF function
 loglsim <- mean(logLam)                     # mean of loglambda
# dse     <- 1.96 * sd(logLam)	              # standard error of loglambda
 dse     <- 1.96 *  apply(logLam, 2, sd)       # to avoid Warning: sd(<matrix>) is deprecated
 CL1     <- c(loglsim - dse, loglsim + dse)  # approx. 95% confidence interval
 lamsim  <- exp(mean(logLam))  		          # simulated stochastic growth rate
 CL2     <- exp(CL1)
## PLOT data
op<-par(mfrow=c(1,2))
hist(logLam, xlab="Log stochastic growth rate",
     col="blue", main=expression(paste("Histogram of log ", lambda[s])) , las=1)
abline(v=CL1, lty=2)
# 
plot(CDFExt,
  type="p", pch=16, col="blue", ylim=c(0,1), las=1,
  main=expression("Extinction CDF"),  
  xlab="Years into the future",
  ylab="Cumulative probability of quasi-extinction")
par(op)
#
vitalsim <- list(
                 detLambda = lam0,
                 stochLambda = c(lambda=lamsim,lc=CL2[1], uc=CL2[2]),
                 logLambdas = c(logLam),
                  CDFExt = CDFExt
                 )
 vitalsim  
} 

