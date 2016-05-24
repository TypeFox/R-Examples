##########################################################################################################
# integrated means and covariences function does univariate and multivariate in one
# function
#
# this function may be used for balanced and unbalanced data so long as the imbalace
# is not extreme
#
# 
# requires:
#
# dat 			a data frame of variables and factors
# data.columns 	 	integer vector indicating which columns in dat contain the measurements
# item.column 	 	integer scaler indicating which column contains the item labels: items 
#			are the top level of variability
#
# returns:
#
# v.within		real p*p matrix where p is the number of variables, estimate of the within
#			fragment covarience matrix
# v.between		real p*p matrix where p is the number of variables - the between items
#			covarience matrix
# n.observations	total number of observations - 1*1 integer
# n.items		number of items - 1*1 integer
# item.n		integer: number of fragments for each item - of length n.fragments
# item.means		real j*p matrix where j is the number of groups and p number of variables
# n.vars		integer: number of contnuous variables
# overall.means		p*1 real vector of means for all variables for all observations
# multivariate		logical T/F indicates whether multivariate
# balanced		logical T/F indicates whether the replication is balanced
# warn.type		character - the type of warning(s) issued - crude for the moment  	
#
# notes:		works by getting the nested means first, then calculating the covariance
#			components from these using weighted estimates
#
# missing values:	crudely handled by mean(x, na.rm=TRUE) to stop any NA's from crashing the
#			function - not very well implemented in this version
##########################################################################################################
two.level.components <- function(dat, data.columns, item.column)
{

# set the warn.type as none - then add in as warnings accrue
warn.type <- "none"

###################### COMMON ############################################################################
# clean the data a bit - get rid of NA rows - crude and may
# lead to cases with n<2 which is tested for later
	if(any(is.na(dat)))
		{
		warning("data contains NAs - cases removed", immediate.=FALSE, call.=FALSE)
		warn.type = "NAs"
		dat <-  dat[complete.cases(dat),]
		}



# definitions and declarations
vars <- names(dat)[data.columns]
n.vars <- length(vars)

multivariate.flag <- TRUE; if(n.vars == 1){multivariate.flag <- FALSE}


n.observations <- nrow(dat)

items <- unique(dat[,item.column])
n.items <- length(items)

item.n <- matrix(0, nrow=n.items, ncol=1); rownames(item.n) <- items; colnames(item.n) <- "n"

item.means <- matrix(0, nrow=n.items, ncol=n.vars); rownames(item.means) <- items; colnames(item.means) <- vars


# means for all cases not affected by imbalance - names attribute akward
# overall mean calculation now even more awkward as comMeans now fussy about
# getting a matrix rather than a vector so colMeans OK if multivariate - otherwise mean
if(multivariate.flag){overall.means <- colMeans(dat[,data.columns], na.rm=TRUE)}
if(!multivariate.flag){overall.means <- mean(dat[,data.columns], na.rm=TRUE);names(overall.means) <- vars}
#overall.means <- colMeans(dat[,data.columns], na.rm=TRUE); if(!multivariate.flag){names(overall.means) <- vars}



# debugging
#assign("kk", overall.means, env=.GlobalEnv)



s.w <- matrix(0, nrow=n.vars, ncol=n.vars); row.names(s.w) <- vars; colnames(s.w) <- vars
s.star <- s.w
##########################################################################################################





###################### multivariate #######################################################################
if(multivariate.flag)
{
	# pick out the observations for each level of the grouping factor then get the mean and the sum of squared deviations (S.w)
	for(ctr in 1:n.items)
		{
		current.item <- dat[dat[,item.column] == items[ctr],]
		item.n[ctr] <- nrow(current.item)
		
		# trap cases with too few replicates to calculate a mean from - fatal if it occurs
		if(item.n[ctr] < 2){stop("too few replicates in an item")}

		item.means[ctr,] <- apply(current.item[,data.columns], 2, mean, na.rm=TRUE)

			
			# for each observation calculate the SS dev component
			for(ctr1 in 1:item.n[ctr])
				{
				s.w <- s.w + (item.n[ctr] * ((as.numeric(current.item[ctr1,data.columns] - item.means[ctr,]) %*% t(as.numeric(current.item[ctr1,data.columns] - item.means[ctr,])))))
				}

		# the between sum of squared deviations between the item means and overall mean - OK for unbalanced at the moment
		s.star <- s.star + (item.n[ctr] * (((item.means[ctr,] - overall.means) %*% t(item.means[ctr,] - overall.means))))
		}
}
##########################################################################################################








###################### univariate ########################################################################
if(!multivariate.flag) # slightly different for univariate case
{

# define the outputs as zero
s.w <- 0
s.star <- 0

	for(ctr in 1:n.items)
		{
		current.item <- dat[dat[,item.column] == items[ctr],]
		item.n[ctr] <- nrow(current.item)



		# differs from multivariate case
		item.means[ctr] <- mean(current.item[,data.columns], na.rm=TRUE)



			# for each observation calculate the SS dev component
			for(ctr1 in 1:item.n[ctr])
				{
				s.w <- s.w + (item.n[ctr] * ((current.item[ctr1, data.columns] - item.means[ctr])^2))
				}
		# the between sum of squared deviations between the item means and overall mean - OK for unbalanced at the moment
		s.star <- s.star + (item.n[ctr] * ((item.means[ctr,] - overall.means)^2))
		}

s.w <- as.matrix(s.w); rownames(s.w) <- vars; colnames(s.w) <- vars
s.star <- as.matrix(s.star); rownames(s.star) <- vars; colnames(s.star) <- vars


}
##########################################################################################################


# both s.w and s.b are pretty well cetain to be correct
# as I have gotten the same values from different independently
# written functions - the final evaluation of C and U is less
# certain - these need a good checking
# give the sums of squared deviations as part of the output
s.w <- s.w * (n.items / n.observations)
s.b <- s.star * (n.items / n.observations)


# set the status of whether the data are balanced between items
if(length(unique(item.n)) > 1){balanced.flag <- FALSE}else{balanced.flag <- TRUE}


U <- (n.items * s.w) / ( n.observations * (n.observations - n.items))
U <- U * (n.observations / n.items) # this bit may be wrong - in so U agrees with previous code

# this may be the correct one
#C <- ((s.star) / (n.observations / n.items * (n.items - 1))) - (s.w / ((n.observations ^ 2 / n.items ^ 2) * (n.observations - n.items)))

# this one may also be wrong - in so U agrees with previous code
#C <- (s.b / (n.items - 1)) - (s.w / ((n.observations^2 / n.items) - n.items))
# this (below) is correct by A&L2004 - thanks to Hanjing Zhang and Colin Aitken for this revision
C <- (s.b / (n.items - 1)) - (s.w / ((n.observations^2 / n.items) - n.observations))

return(new("compcovar", v.within=U, v.between=C, n.observations=n.observations, n.items=n.items, item.n=item.n, item.means=item.means, n.vars=n.vars, overall.means=overall.means, multivariate=multivariate.flag, balanced=balanced.flag, s.within=s.w, s.between=s.b, warn.type=warn.type))
}

