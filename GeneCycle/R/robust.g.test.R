### robust.g.test.R  (2008-07-02)
###
###    Robust detection of periodic time series
###
### Copyright 2005-8, 2008-7 Miika Ahdesmaki
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
# robust.g.test.R
######################################################################
######################################################################

# public functions

#robust.g.test <- function(y, index, perm = FALSE, x, noOfPermutations = 5000) 
robust.g.test <- function(y, index, perm = FALSE, x, noOfPermutations = 300, algorithm=c("rank", "regression"), t) 
{
  # ASSUMPTION: Each time series is of the same length.
  algorithm = match.arg(algorithm)

  # 'y' is the matrix consisting of the spectral estimates
  # as column vectors. 
  
  y <- as.matrix(y)
  
  # 'index' is an index to the spectral
  # estimates that is to be used in the testing for periodicity.
  # If 'index' is missing, the maximum component of the spectral
  # estimate is used in testing (regardless of the frequency 
  # of this maximum). 
  # NOTE: 'index' is unrelevant for the regression approach; only 
  # the maximum component will be tested in this function. For testing
  # a given frequency in the regression approach, please see the
  # file robust.spectrum.R
  if (algorithm == "regression" && !missing(index)) {rm(index)}
  
  # 'perm' is a boolean. If 'perm' is FALSE, 
  # a simulated distribution for the g-statistic is used. 
  # If per 'perm' is TRUE, permutation  tests are used to find 
  # the distribution of the g-statistic for each time series
  # separately. 
  # NOTE: permutations will be used for the regression approach
  # regardless of the given state of this variable.
  if (algorithm == "regression") {perm = TRUE}
  
  #'x' must include the original time series for 
  # the permutation alternative (not needed if permutation 
  # tests are not used).
  
  # 'noOfPermutations' tells the number of permutations that 
  # are used for each time series (default = 300).

  # Check for the regression based approach if 'x' exists
  # and whether the time series are of even or odd length
  if (algorithm == "regression") {
    if(missing(x) || missing(t)){stop("Original time series x and sampling time vector t must be specified")}
    if(nrow(x) %% 2 == 0){
      isEven = TRUE
    }else{isEven = FALSE}
    gstats <- g.statistic(y, isEven=isEven)
  }else{
  gstats <- g.statistic(y, index)	# Evaluate the g-statistic
  }						# for each spectrum

  ######################################################
  # First we consider the approach where we simulate the
  # distribution of the g-statistic under the null hypothesis
  ###########################################################
  if(perm == FALSE){
  	tSLength <- nrow(y)	# length of one time series
  	
	if (missing(index)){
  		filename <- paste("g_pop_length_",tSLength,".txt",sep="")
  	}
  	else{
  		filename <- paste("g_pop_length_",tSLength,"indexed.txt",sep="")
  	}
  	# If a null hypothesis distribution does not already exist
  	# for this time series length, generate a new one.
  	if (file.exists(filename)==FALSE){
  		gPopCreate(tSLength, filename, index)
  	}
  	temp <- as.matrix(read.table(filename))
  	gPop <- matrix(temp,length(temp),1)
  	
	
	gPop.dens <- density(gPop)	# Kernel density estimate
  	# gPop.dens$x is now the g-statistic vector and gPop.dens$y
  	# the estimated distribution value-vector.
  	# N.B. that gPop.dens$y is now evenly spaced
  	# (..$x[n]-..$x[n-1]=constant), making computations easier. 
  	# sum(gPop.dens$y)*(gPop.dens$x[2] - gPop.dens$x[1]) = 1
  	# (approximately)
  	
  	# Now let us find the p-values:
  	pvals <- 0
  	for (k in 1:length(gstats)){
  		pvals[k] <- sum(gPop.dens$y[gPop.dens$x > gstats[k]])
  	}
  	# multiplier <- (gPop.dens$x[11] - gPop.dens$x[1])/10
  	multiplier <- (gPop.dens$x[2] - gPop.dens$x[1])
  	pvals <- pvals * multiplier
  	return(pvals)
  }
  #########################################################
  # Then we consider the permutation test based alternative
  # for finding p-values
  # REMEMBER: This is a very slow and computing intensive
  # approach!
  ###########
  else{
  	if ( missing(x) ) 
	  stop("Original time series x must be specified")

	tSLength <- nrow(y)	# length of one time series
  	tSNumber <- ncol(y)	# Number of time series
  	pvals <- 0			# Initialise p-values
  	gs <- 0			# Initialise g-statistics
  	for (tS in 1:tSNumber){
  		temp <- matrix(0,length(x[,tS]),noOfPermutations)
  		for (perm in 1:noOfPermutations){
  			# Resample the time series...
  			temp[,perm] <- myresample(x[,tS])
  		}
	   if(algorithm=="rank"){
  		# ...And evaluate the Pearson's spectral estimator
  		# for each permutation
  		specTemp <- robust.spectrum(temp) 
  		# Evaluate g-statistics for each spectral estimate
  		gs <- g.statistic(specTemp, index)
           }else{
		# evaluate the regression based estimator
  		specTemp <- robust.spectrum(x = temp, algorithm = "regression", t=t) 
  		gs <- g.statistic(y = specTemp, isEven = isEven)		
           }

  		# Form an estimate of the distribution of the 
  		# g-statistic
  		gsDens <- density(gs)
  
  		# Finally, estimate the p-values
  		pvals[tS] <- sum(gsDens$y[gsDens$x > gstats[tS]])
  		multiplier <- (gsDens$x[2] - gsDens$x[1])
  		pvals[tS] <- pvals[tS] * multiplier
  	}
  	
  	return(pvals)	
  }	# If..else ends
}	# Function ends



##########################################################
###Returns the g-statistics of spectra####################
##########################################################
g.statistic <- function(y, index, isEven=FALSE)
{
  if(missing(index)){
	nrowsy <- nrow(y)
	if(isEven){
  	maxVals <- apply(y[c(-1,-nrowsy), ,drop=FALSE],2,max)	# Maximum component
  	sumVals <- apply(y[c(-1,-nrowsy), ,drop=FALSE],2,sum)	# Sum of components
  	return(maxVals / sumVals)	# Return the ratios
	}else{
  	maxVals <- apply(y[-1, ,drop=FALSE],2,max)	# Maximum component
  	sumVals <- apply(y[-1, ,drop=FALSE],2,sum)	# Sum of components
  	return(maxVals / sumVals)	# Return the ratios
	}
  }
  else{
  	vals <- y[index,]	# Picks the index:th component of each spectrum
  	sumVals <- apply(y[-1,],2,sum)
  	return(vals / sumVals)
  }

# note: y[-1,] has to be used instead of y to exclude zero frequency
# The pi-frequency is never included in the rank-spectra (due to
# the implemented interpolation and exclusion of the pi-frequency)
# The pi-frequency is included in the regression spectra (when N 
# even) but is not considered in calculating the g-statistic.
# changed 3 July 2008
}


######################################################################

# private functions


##########################################################
####Inner function : Create a g-statistic "population"####
##########################################################
gPopCreate <- function(tSLength, filename, index, size = 10000)
{
  # This function generates a file that includes the 
  # g-statistics of random time series to estimate the null 
  # hypothesis distribution of g.
  temp <- rnorm(tSLength*size)
  # print("Creating...")
  # Create a matrix of "size" time series, time series 
  # length "tSLength"
  temp <- matrix(temp,tSLength,size) 	
  spectra <- robust.spectrum(temp)
  gStats <- g.statistic(spectra, index)# Evaluate the g-statistic
  					# for each spectrum
  write(gStats, filename)
}

##########################################################
####Inner function : Resampling for permutation tests ####
##########################################################
myresample <- function(x, size, ...)
{
  if(length(x) <= 1) { 
  	if(!missing(size) && size == 0) x[FALSE] 
  	else x
  } 
  else sample(x, size, ...)
}

######################################################################
