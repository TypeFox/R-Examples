getNormalLHS = function(N, mean, sd, aboveZeroOnly = TRUE, belowOneOnly = FALSE) {
  # returns a LHS from the normal distibution with
  # mean and sd: splits the normal distribution
  # into N equiprobable sections and then returns a
  # shuffled vector containing a sample from each section

  # N random numbers between 0 and 1/N
  randProb = runif(N)/N
  #print(randProb)

  # lower probability bounds for N equally weighted
  # probability sections
  sectionLowerBounds = 0:(N - 1)/N
  #print(sectionLowerBounds)

  # generates random probability taken from each
  # of the N equally weighted probability sections
  rand = randProb + sectionLowerBounds
  #print(rand)

  # reverse quantile function to get random samples
  # from each equally weighted section of the normal
  # probability distribution with mean and sd
  LHsample = qnorm(rand, mean, sd)
  #print("preRemoveNegatives")
  #print(LHsample)

  # if any values are below 0, set to 0
  ######## WARNING!!! : YOU MAY NOT WANT THIS!!! #############
  if (aboveZeroOnly) {
    LHsample[LHsample < 0] <- 0
  }

  # if any values are above 1, set to 1
  ######## WARNING!!! : YOU MAY NOT WANT THIS!!! #############
  if(belowOneOnly){
    LHsample[LHsample > 1] <- 1
  }

  #print("postRemoveNegatives")
  #print(LHsample)

  # takes a sample from LHsample without replacement
  shuffledLHS = sample(LHsample)

  # return shuffled latin hypercube sample
  shuffledLHS
}

getUniformLHS = function(N, min, max) {
  # returns a LHS from the uniform distibution with
  # min and max: splits the uniform distribution
  # into N equiprobable sections and then returns a
  # shuffled vector containing a sample from each section

  # N random numbers between 0 and 1/N
  randProb = runif(N)/N
  #print(randProb)

  # lower probability bounds for N equally weighted
  # probability sections
  sectionLowerBounds = 0:(N - 1)/N
  #print(sectionLowerBounds)

  # generates random probability taken from each
  # of the N equally weighted probability sections
  rand = randProb + sectionLowerBounds
  #print(rand)

  # reverse quantile function to get random samples
  # from each equally weighted section of the uniform
  # probability distribution with min and max
  LHsample = qunif(rand, min, max)
  #print(LHsample)

  # takes a sample from LHsample without replacement
  shuffledLHS = sample(LHsample)

  # return shuffled latin hypercube sample
  shuffledLHS
}

## Summary SE Summarizes data and is Taken from: http://www.cookbook-r.com/Manipulating_data/Summarizing_data/

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE, conf.interval = 0.95, .drop = TRUE) {

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function(x, na.rm = FALSE) {
    if (na.rm)
      sum(!is.na(x))
    else length(x)
  }

  # This is does the summary; it's not easy to understand...
  datac <- plyr::ddply(data, groupvars, .drop = .drop,
                       .fun = function(xx, col, na.rm) {
                         c(N = length2(xx[, col], na.rm = na.rm),
                           mean = mean(xx[, col], na.rm = na.rm),
                           sd = sd(xx[, col], na.rm = na.rm),
                           max = max(xx[, col], na.rm = na.rm),
                           min = min(xx[, col], na.rm = na.rm),
                           q1 = quantile(xx[, col], na.rm = na.rm, 1/4),
                           q3 = quantile(xx[, col], na.rm = na.rm, 3/4),
                           median =  median(xx[, col], na.rm = na.rm))
                       },
                       measurevar, na.rm)

  # Rename the "mean" column
  datac <- plyr::rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd/sqrt(datac$N) # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + 0.5, datac$N - 1)
  datac$ci <- datac$se * ciMult

  colnames(datac)[7] <- "q1"
  colnames(datac)[8] <- "q3"

  return(datac)
}

xLimsBarPlot = function(numBars, spacing = 0.2, padding = 0.04){
  xhigh = numBars*(1+spacing)*(1+padding)
  diff = xhigh - numBars*(1+spacing)
  xlow = (0+spacing)-diff
  #print(c(xlow, xhigh))
  return(c(xlow, xhigh))
}

### Model ###
model = function(sales, fakePercent, deathRate) {
  #fakeProportion = fakePercent/100
  fakeProportion = fakePercent
  deaths = sales * fakeProportion * deathRate

  # add sum to end of vector
  deaths = c(deaths, sum(deaths))

  # return deaths
  deaths
}

# returns the standard deviation of a uniform distribution
# with a given min and max
stdUniform = function(min, max){
  std <- sqrt((1/12)*(max-min)^2)
}

# counterfeitPRCC = function(x, sort.results = FALSE, sort.abs = FALSE) {
# 	# This code was adopted from Eili Klein <klein@cddep.org>
#
# 	# N = number of simulations
# 	N = length(x[, 1])
# 	#print("N:")
# 	#print(N)
#
# 	# k = number of input parameters + 1 output var
# 	k = length(x[1, ])
# 	#print("K:")
# 	#print(k)
#
# 	# Rank each input parameter
# 	r = matrix(NA, nrow = N, ncol = k)
# 	for (i in 1:k) {
# 		#r[, i] = rank(ties.method = "min",x[, i])
# 		r[, i] = rank(x[, i])
# 	}
#
# 	# save output ranks
# 	outputvar = r[, k]
#
# 	# r is the matrix where each column (1 for each input parameter)
# 	# contains ranks of the simulation values for that parameter
# 	r = r[, -k]
# 	k = k - 1
#
# 	# If two of the input parameters have exactly the same ranking for
# 	# every run, then only one of the parameters should be used in the
# 	# calculation of PRCC
# 	dropcols = 0
# 	for (i in 1:(k - 1)) {
# 		dups = seq(1, k)
# 		for (j in 1:i) dups[j] = 0
# 		for (j in (i + 1):k) {
# 			a = which(r[, i] == r[, j])
# 			if (length(a) != N) {
# 				dups[j] = 0
# 			}
# 		}
#
# 		# keep track of duplicate columns
# 		for (j in 1:k) {
# 			if (dups[j] > 0) {
# 				dropcols = c(dropcols, j)
# 			}
# 		}
# 	}
#
# 	# remove 0 as a dropcol
# 	dropcols = dropcols[-1]
#
# 	# get unique column numbers
# 	dropcols = unique(dropcols)
#
# 	# sort dropcol list
# 	dropcols = sort(dropcols)
#
# 	# reverse list of columns to be dropped so that you
# 	# can drop them using the # as an index and it wont
# 	# interfere with subsequent drops
# 	dropcols = rev(dropcols)
#
# 	# drop duplicate columns
# 	if (length(dropcols) > 0) {
# 		for (i in 1:length(dropcols)) {
# 			r = r[, -dropcols[i]]
# 		}
# 	}
#
# 	# adjust K if you droped any columns
# 	k = length(r[1, ])
#
# 	# bind output rankings to input rankings matrix
# 	r = cbind(r, outputvar)
#
# 	### Generate C Matrix ###
# 	mu = (1 + N)/2
# 	C.ij = matrix(NA, nrow = k + 1, ncol = k + 1)
# 	for (i in 1:(k + 1)) {
# 		for (j in 1:(k + 1)) {
# 			C.ij[i, j] = sum((r[, i] - mu) * (r[, j] - mu))/sqrt(sum((r[, i] -
# 				mu)^2) * sum((r[, j] - mu)^2))
# 		}
# 	}
#
# 	B = rms::matinv(C.ij)
#
# 	gamma.ij = rep(0, k)
# 	t.iy = rep(0, k)
# 	p.iy = rep(0, k)
#
# 	for (i in 1:(k)) {
# 		# the PRCC between the ith input parameter and the outcome variable
# 		gamma.ij[i] = -B[i, k + 1]/sqrt(B[i, i] * B[k + 1, k + 1])
#
# 		# the significance of a nonzero PRCC is tested by computing t.iy
# 		# the distribution of t.iy approximates a students T with N-2 degrees of freedom
# 		t.iy[i] = gamma.ij[i] * sqrt((N - 2)/(1 - gamma.ij[i]))
#
# 		p.iy[i] = 2 * pt(-abs(t.iy[i]), df = N - 2)
# 	}
#
# 	# account for dropped columns
# 	dropcols = sort(dropcols)
# 	if (length(dropcols) > 0) {
# 		for (i in 1:length(dropcols)) {
# 			if (dropcols[i] > length(gamma.ij)) {
# 				# append to end
# 				gamma.ij = c(gamma.ij[1:(dropcols[i] - 1)], 0)
# 				t.iy = c(t.iy[1:(dropcols[i] - 1)], "dropped")
# 				p.iy = c(p.iy[1:(dropcols[i] - 1)], "-")
# 			} else {
# 				# insert in location
# 				gamma.ij = c(gamma.ij[1:(dropcols[i] - 1)], 0, gamma.ij[dropcols[i]:(length(gamma.ij))])
# 				t.iy = c(t.iy[1:(dropcols[i] - 1)], "dropped", t.iy[(dropcols[i]):(length(t.iy))])
# 				p.iy = c(p.iy[1:(dropcols[i] - 1)], "-", p.iy[(dropcols[i]):(length(p.iy))])
# 			}
# 		}
# 	}
#
# 	# calculate absolute value column to be used for sorting
# 	gamma.ij.abs = abs(gamma.ij)
#
# 	# create output dataframe
# 	vals = data.frame(cbind(gamma.ij, t.iy, p.iy, gamma.ij.abs), row.names = names(x[1:(length(x[1,
# 		]) - 1)]))
#
# 	if (sort.results == TRUE) {
# 		if (sort.abs == TRUE) {
# 			vals = vals[order(vals$gamma.ij.abs, decreasing = TRUE), ]
# 		} else {
# 			vals = vals[order(vals$gamma.ij, decreasing = TRUE), ]
# 		}
# 	}
#
# 	# drop absolute value column
# 	vals = vals[, -4]
#
# 	vals
#
# }

