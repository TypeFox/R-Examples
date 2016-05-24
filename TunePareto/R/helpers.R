#################################################################
# Helper functions
#################################################################


# Iteratively calculate all combinations of variable assignments,
# where <numValues> is a vector specifying the number of choices
# for each of the variables.
# <currentCombination> is the last used combination from which
# the next new combination is calculated.
#
# Returns the next combination or NULL if all combinations have
# been tested.
nextCombination <- function(numValues, currentCombination)
{	
	while(TRUE)
	{
		for(i in (length(currentCombination)):1)
		{
			currentCombination[i] <- currentCombination[i] + 1	
			if(currentCombination[i] <= numValues[i])
			{
				return(currentCombination)	
			}
			else
			{
				if (i==1)
					return(NULL)	
				else
					currentCombination[i] <- 1
			}		
		}		
	}
}

# Retrieves a data frame of possible parameter combinations
# based on the above method.
# <parameterRanges> is a list of lists, where each of the
# sub-lists specifies the possible values a parameter can take.
allCombinations <- function(parameterRanges)
{
  combination <- rep(1, length(parameterRanges))
  
  # build a data frame from the first combination
  res <- list(mapply(function(parameters, index)
                             {
                              parameters[[index]]
                             },
                             parameterRanges, combination,
                             SIMPLIFY=FALSE))

  numValues <- sapply(parameterRanges,length)

  j <- 1
  while (TRUE)
  {
    j <- j + 1
    combination <- nextCombination(numValues, combination)
    if (is.null(combination))
    # no more combinations
      break
      
    # append new combination to data frame
    res[[j]] <- mapply(function(parameters, index)
                             {
                              parameters[[index]]
                             },
                             parameterRanges, combination,
                             SIMPLIFY=FALSE)
   }
   return(res)
}

# Draws a sample of <N> combinations from <parameterRanges>
# using the specified sampling method.
sampleCombinations <- function(parameterRanges, N, method=c("uniform","halton","niederreiter","sobol"))
{
  # first, draw a set of values in [0,1]
  numbers <- switch(match.arg(method,c("uniform","halton","niederreiter","sobol")),
                    uniform = runif(n=length(parameterRanges)*N),
                    halton = halton.sample(n=N, dimension=length(parameterRanges)),
                    niederreiter = {
                                      require(gsl)
                                      pt <- qrng_alloc(type="niederreiter_2",dim=length(parameterRanges))
                                      
                                      # random initialization
                                      qrng_get(pt,sample(1:1000,size=1))
                                      t(qrng_get(pt, N))
                                      
                                   },
                    sobol =        {
                                      require(gsl)
                                      pt <- qrng_alloc(type="sobol",dim=length(parameterRanges))
                                      
                                      # random initialization
                                      qrng_get(pt,sample(1:1000,size=1))
                                      t(qrng_get(pt, N))
                                   }                                   
                   )
  numbers <- matrix(numbers,ncol=N)
  
  # apply the [0,1] values to the ranges of the parameters
  res <- apply(numbers,2,function(col)
                {
                  x <- mapply(function(val,range)
                  {
                    if (is.interval(range))
                    # interval => scaling
                    {
                     range$lower + (range$upper-range$lower) * val
                    }
                    else
                    # discrete parameter => round to "grid"
                    {
                      if (val == 0)
                        val <- 0.01
                      range[ceiling(val * length(range))]
                    }
                  }, col, parameterRanges, SIMPLIFY=FALSE)
                  names(x) <- names(parameterRanges)
                  x
                })
  return(res)                
}

# Latin Hypercube sampling of <N> combinations from the
# ranges supplied in <parameterRanges>
latinHypercube <- function(parameterRanges, N)
{
  vals <- data.frame(lapply(parameterRanges,function(range)
          {
            if (is.interval(range))
            # partition interval into <N> parts and take
            # one value from each part
            {
              intSz <- (range$upper - range$lower)/N
              lower <- cumsum(c(range$lower,rep(intSz,N-1)))
              return(sample(sapply(lower,function(x)runif(min=x,max=x+intSz,n=1)),
                            size=N,replace=FALSE))
            }
            else
            # ensure that each discrete value occurs equally often (+-1)
            {
              if (N %% length(range) == 0)
                reps <- N %/% length(range)
              else
                reps <- N %/% length(range) + 1
              return(sample(rep(range,reps),size=N,replace=FALSE))
            }
          }),stringsAsFactors=FALSE)
  return(lapply(1:nrow(vals),function(i)
            {
              x <- as.list(vals[i,])
              names(x) <- names(parameterRanges)
              x
            }))
}

# Creates an interval object from
# a lower and an upper bound.
as.interval <- function(lower,upper)
{
  l <- list(lower=lower,upper=upper)
  class(l) <- "Interval"
  return(l)
}

# Check whether an object is an interval.
is.interval <- function(x)
{
  return(inherits(x,"Interval"))
}

