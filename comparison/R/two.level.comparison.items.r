##########################################################################################################
# 
# Makes a comparison.items object. A comparison items object has the data from the item to be 
# compared - typically for any specific comparison problem there will be two of these objects
# 
# requires:
#
# dat 			a data frame of variables and factors
# data.columns 	 	integer vector indicating which columns in dat contain the measurements
#
# returns:
#
# item.means		real p array where p is the number of variables
# item.n		integer: number of replicate observations for each item
# n.vars		integer: number of contnuous variables
# multivariate		logical T/F indicates whether multivariate or not
# observed		the raw observations in the form of an i*p matrix where there are i replicated
#			observations and p variables
#
# missing values:	traps and removes NAs - issues warnings about NAs
##########################################################################################################
two.level.comparison.items <- function(dat, data.columns)
{

# set the warn.type as none - then add in as warnings accrue
warn.type <- "none"

# clean the data a bit - get rid of NA rows - crude and may
# lead to cases with n<2 which is tested for later
	if(any(is.na(dat)))
		{
		warning("data contains NAs - cases removed", immediate.=FALSE, call.=FALSE)
		warn.type = "NAs"
		dat <-  dat[complete.cases(dat),]
		}


# requires a cast to matrix for some reason
observed <- as.matrix(dat[, data.columns], rownames.force=TRUE)

# simple counts of the data
item.means <- apply(observed, 2, mean)
n.replicates <- nrow(observed)
n.vars <- length(data.columns)

# test for, and set the multivariate flag
if(n.vars > 1){multivariate.flag <- TRUE} else{multivariate.flag <- FALSE}

# tidy up the rownames property of the observations
rownames(observed) <- as.character(1:n.replicates)


return(new("compitem", item.means=item.means, n.replicates=n.replicates, n.vars=n.vars, multivariate=multivariate.flag, observed=observed, warn.type=warn.type))
}

