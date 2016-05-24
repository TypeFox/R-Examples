### robust.spectrum.R  (2009-04-09)
###
###    Robust estimates of spectra of time series
###
### Copyright 2005, 2009 Miika Ahdesmaki
###
###
###
### This file is part of the `GeneCycle' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA




######################################################################
######################################################################
# robust.spectrum.R
######################################################################
######################################################################

# public functions


#################################################################
############# The spectra of multiple time series ###############
#################################################################
robust.spectrum <- function(x,algorithm=c("rank", "regression"), t, periodicity.time=FALSE, noOfPermutations=300) 
# 'x' is the matrix of time series as column vectors
# 
# For the regression approach, the time vector 't' of sampling times 
# must be given.
#
# The variable 'frequency' is by default FALSE and will only affect the
# regression based approach. If it has a numerical (positive) value, the 
# given frequency will be used in testing for periodicity (only for the
# regression based method)
# 'noOfPermutations' is the number of permutations used in the permutation tests - 
# however, it only affects the regression based method when periodicity.time
# is given some numerical value
#
#
# The output of the function will be the estimated spectra (to be 
# further processed by robust.g.test) or, for the regression 
# based method _and_ known periodicity frequency, the p-values.
{
  algorithm = match.arg(algorithm)

  if (algorithm == "rank"){
    ## BEGIN RANK##
    # x should be a matrix consisting of the time series as column
    # vectors.
    if (is.longitudinal(x))
    {
      if (has.repeated.measurements(x) )
        stop("Multiple measurements per time point are not supported")
    
      if (!is.equally.spaced(x))
        warning("Function assumes equally spaced time points")
    }
    else
    {
      x <- as.matrix(x)
    }

    noTS <- ncol(x) #number of time series
    y <- matrix(NA,nrow(x),ncol(x))
    for (h in 1:noTS)
    {
      temp <- x[,h]
      y[,h] <- robust.spectrum.single(temp)
    }  
    return(y)
  ## END RANK##


  }else{
  ## BEGIN REGRESSION APPROACH
    if(missing(t)){
      stop("The sampling times must be given for the regression based approach.")}
    if (is.longitudinal(x))
    {
      if (has.repeated.measurements(x) )
        stop("Multiple measurements per time point are not supported")
    
      if (!is.equally.spaced(x))
        warning("Please give the time indices in a separate vector t!")
    }
    else
    {
      x <- as.matrix(x)
    }

	require(MASS)
	# library(lqs)

    # two cases: search all frequencies (after which apply robust.g.test)
    # or just one given frequency (for which p-values will be evaluated here)
    if(periodicity.time){
    ############ KNOWN FREQUENCY, REGRESSION BASED APPROACH
	noOfGenes = ncol(x)
	Nt = nrow(x)
	if(!(Nt==length(t))){stop("Length of t must be equal to the number of rows in x")}
	timeVector = (t - t[1]) / (t[Nt] - t[1]) * (Nt-1) # "index"-vector with non-integer values
	Fs = (timeVector[2]-timeVector[1]) / (t[2]-t[1]) # 1/The mean sampling time
	f = (1/Fs)/periodicity.time # the normalised frequency at which periodicity will be sought
	sine = sin(2*pi*f*timeVector)   # The sine we use in regression
	cosine = cos(2*pi*f*timeVector) # The cosine we use in regression

	A = matrix(0, 2,noOfGenes)

	for (gene in 1:noOfGenes){
		temp = rlm(x[,gene] ~ sine + cosine, method='MM')
    		A[,gene] = matrix(temp$coefficients[2:3],2,1)
	}
	A=A[1,]^2 + A[2,]^2
	
	distri = matrix(0,noOfPermutations,noOfGenes)
	for (kierros in 1 : noOfPermutations){
	    permu = sample(Nt)
	    for (gene in 1:noOfGenes){
	        sarja = x[permu,gene]
	        temp = rlm(sarja ~ sine + cosine, method='MM')
	        distri[kierros,gene] = temp$coefficients[2]^2 + temp$coefficients[3]^2
	    }
	}

	density_size = 512	#density()-default size
	xDens = matrix(0, density_size,noOfGenes)
	yDens = matrix(0, density_size,noOfGenes)
	for (gene in 1 : noOfGenes){
	    #[FF_rob(:,gene),Xi_rob(:,gene)] = density(distri[,gene])
	    temp = density(distri[,gene], from=0)
	    xDens[,gene] = temp$x
	    yDens[,gene] = temp$y
	}
	pval = matrix(0,noOfGenes,1)
	for (gene in 1:noOfGenes){ #calculate p-values for each test/gene
		pval[gene] = sum(yDens[A[gene] < xDens[, gene], gene])
		pval[gene] = pval[gene] * (xDens[2,gene]-xDens[1,gene]) #assumption:...
	    #...density() returns equally sampled densities
	}
	pval[pval>1] = 1	# Just to be sure
	print("Please note, returning p-values!")
	return(pval)
    }else{
    ############## UNKNOWN FREQUENCY, REGRESSION BASED APPROACH
	noOfGenes = ncol(x)
	Nt = nrow(x)
	if(!(Nt==length(t))){stop("Length of t must be equal to the number of rows in x")}
	timeVector = (t - t[1]) / (t[Nt] - t[1]) * (Nt-1) # "index"-vector with non-integer values
	Fs = (timeVector[2]-timeVector[1]) / (t[2]-t[1]) # 1/The mean sampling time
	if(Nt %% 2 == 0){
      		isEven = TRUE
	}else{isEven = FALSE}
	k = 0 : floor((Nt)/2)		# if Nt is even, the component Nt/2 will be considered separately
	K = length(k)			# 
	kronpro = matrix(k %x% timeVector, Nt, K)
	sines = sin(2*pi*kronpro/Nt)	# The sines we use in regression
	cosines = cos(2*pi*kronpro/Nt)	# The cosines we use in regression
	
	# the coefficients/wieghts for the sines and cosines will be collected in A
	#
	sinCoeffs=matrix(0,K,noOfGenes)
	cosCoeffs = matrix(0,K,noOfGenes)

	# Next we first fit the sines (and cosines) to the signals without residual subtraction (to get an initial estimate). Then, based on the initial estimate, we again fit the sines in decreasing order of power of the initial fit (i.e. the frequency, which has the most power in the initial estimate, is fitted first and so on) and always fit only to the residuals of the previous round (the first frequency is fit to the original signal).
	for (gene in 1:noOfGenes){
		temp = fit.DC(x[,gene])
		cosCoeffs[1,gene] = temp$coeffs[1]	# intercept/DC-level
		sinCoeffs[1,gene] = temp$coeffs[2]	# 0
		if(isEven){
			for (frek in 2:(K-1)){
				temp = fit.freqs(x[,gene],sines[,frek],cosines[,frek])
				cosCoeffs[frek,gene] = temp$coeffs[1]
				sinCoeffs[frek,gene] = temp$coeffs[2]
			}
			temp = fit.pi(x[,gene],cosines[,K])
			cosCoeffs[K,gene] = temp$coeffs[1]	# pi-frequency
			sinCoeffs[K,gene] = temp$coeffs[2]	# 0
			
		}else{
			for (frek in 2:(K)){
				temp = fit.freqs(x[,gene],sines[,frek],cosines[,frek])
				cosCoeffs[frek,gene] = temp$coeffs[1]
				sinCoeffs[frek,gene] = temp$coeffs[2]
			}
		}
	}

	# ... and then fit to residuals
	tempSpectrum = cosCoeffs^2 + sinCoeffs^2
	#print(tempSpectrum)
	sinCoeffs = matrix(0,K,noOfGenes)
	cosCoeffs = matrix(0,K,noOfGenes)	
	# B = list(sinCoeffs=matrix(0,K,noOfGenes), cosCoeffs = matrix(0,K,noOfGenes))
	for (gene in 1:noOfGenes){
		sortSpec = sort(tempSpectrum[,gene], decreasing = TRUE, index.return = TRUE)
		temp = fit.DC(x[,gene])
		cosCoeffs[1,gene] = temp$coeffs[1]	# intercept/DC-level
		sinCoeffs[1,gene] = temp$coeffs[2]	# 0
		RESID = temp$resid
		if(isEven){
			for (frek in sortSpec$ix){
				if(frek==1 || frek == K){next}
				temp = fit.freqs(RESID,sines[,frek],cosines[,frek])
				cosCoeffs[frek,gene] = temp$coeffs[1]
				sinCoeffs[frek,gene] = temp$coeffs[2]
				RESID = temp$resid
			}
			temp = fit.pi(RESID,cosines[,K])
			cosCoeffs[K,gene] = temp$coeffs[1]	# pi-frequency
			sinCoeffs[K,gene] = temp$coeffs[2]	# 0
			
		}else{
			for (frek in sortSpec$ix){
				if(frek==1){next}
				temp = fit.freqs(RESID,sines[,frek],cosines[,frek])
				cosCoeffs[frek,gene] = temp$coeffs[1]
				sinCoeffs[frek,gene] = temp$coeffs[2]
				RESID = temp$resid
			}
		}	
		
		
	}
	y = cosCoeffs^2 + sinCoeffs^2
	#print(y)
    	return(y)
    }

  ############## END REGRESSION APPROACH
  }

}

######################################################################

# private functions


#################################################################
# This function implements the robust, rank-based spectral#######
# estimator introduced in Pearson et al. 2003####################
#################################################################
robust.spectrum.single <- function(x) 
{
  
  if( is.constant.single(x) ) warning("Constant time series!")
  
  ##############################################
  # Some adjustable parameters
  ##############################################
  # Length of the zero-padded one-sided "rho"(=Rsm)
  zp <- 2*length(x)
  
  # Length of the original sequence
  n <- length(x)
  
  # Let us define the maximum lag for the correlation coefficient:
  maxM <- n-2
  
  
  ##############################################
  # Correlation coefficient
  ##############################################
  # Reserve space
  Rsm <- matrix(NA, nrow = maxM+1, ncol = 1)


  missing <- is.na(x)
  nonMissing <- !missing
  mi <- which(missing)     # indices for missing values
  nmi <- which(nonMissing) # indices for nonmissing values
  
  # Mean removal
  x[nmi] <- x[nmi] - mean(x[nmi])
  
  # Run through all the lags
  for (lags in 0:maxM)
  {
    # Modified Spearman's method
    indexes <- 1:(n-lags)	# Initial indices
    ends <- length(nonMissing)
  
    # Values in both the original and shifted vectors must be present:
    temp <- (nonMissing[1:(ends-lags)] + nonMissing[(lags+1):ends]) >= 2
    indexPresent <- which(temp)
    indexes <- indexes[indexPresent]	# The indices that are present in
  							# both sequences
    Rsm[lags+1] <- ifelse(
                     length(indexes)<=1 , 
		     0 ,
		     spearman(x[indexes],x[indexes+lags],n))
  }
  
  # Zero-padding
  Rsm[(length(Rsm)+1):zp] <- 0
  fftemp <- fft(Rsm)
  
  # The following implementation is as in (Ahdesmäki, Lähdesmäki et al., 2005)
  Ssm <- abs( 2*Re(fftemp) - Rsm[1] )
  Ssm <- Ssm[1:floor(length(Ssm)/2)]
  
  # Return the spectral content, frequencies [0,pi)
  return(Ssm)
}


#################################################################
########Inner function:Spearman's correlation coefficient########
#################################################################


spearman <- function(x, y, N, version=c("builtin", "miika") )
{
  #cat (paste("DEBUG: ", length(x), length(y), N, "\n") )
  
  version <- match.arg(version)
  
  if (version == "builtin")
  {
      rho <- cor(x, y, method="spearman" ) * length(x)/N
  }
  
  if (version == "miika")
  {
      Km <- length(x)
      sx <- sort(x); ix <- order(x);
      sy <- sort(y); iy <- order(y);

      rx <- matrix(0, 1, Km)
      ry <- rx
      rx[ix] <- (1:Km)
      ry[iy] <- (1:Km)

      rho <- 12/(N*(Km^2-1)) * (rx-(Km+1)/2) %*% (t(ry) - (Km+1)/2)
  }
  
  return(rho)
}

#############################################
# check that the two versions are the same
#
# x1 <- rnorm(10)
# x2 <- rnorm(10)
#
# spearman(x1, x2, 111)
# spearman(x1, x2, 111, version="miika")
#############################################




#################################################################
#################################################################
######## Private functions for the regression approach   ########
#################################################################
#################################################################

#################################################################
########Private function: Regression approach, DC-level  ########
#################################################################
fit.DC <- function(x){
	temp = rlm(x ~ 1, method='MM')
	B = list(coeffs = c(2 * temp$coefficients[1],0), resid = temp$residuals)
	return(B)
}
#################################################################
########Private function: Regression approach, Fourier freqs ####
#################################################################
fit.freqs <- function(x,sine,cosine){
	temp = rlm(x ~ sine + cosine, method='MM')
	B = list(coeffs = temp$coefficients[2:3], resid = temp$residuals)
	return(B)
}
#################################################################
########Private function: Regression approach, pi-frequency #####
#################################################################
# NOTE: this function should only be called if the time series length is even
fit.pi <- function(x,cosine){
	# note: for pathological sampling (eg. 0, 1.5, 2.5, 3.5, 4.5,...) this implementation is bad, since we are not fitting the sinusoidal (for integer, i.e. ordinary, sampling the sinusoidal term is always zero)
	temp = rlm(x ~ cosine, method='MM')
	B = list(coeffs = c(2 * temp$coefficients[2],0), resid = temp$residuals)
	return(B)
}
#################################################################
