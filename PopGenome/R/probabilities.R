###############################################################################
#
# FUNCTION: calc_probabilities, calc_quantile
#
# This function compares the observed value with the data which was
# gained by coalescent simulation and calculates probabilties for the occurence
# of the observed values under the hypothesis the coalescent simulation was run
# 
#
# FUNCTION CALLS: 	
#					
#
# PARAMETERS:
#		cloci:		current loci
#
#		niter:		number of iterations for current loci
#
#		csStats:	a matrix containing summary statistics for current loci
#
#		obsVal:		a matrix with observed value, given by the user
#
#		lData:		an object of class locStats, for saving all computed data
#
#
# RETURN VALUES:	an object of class locstats, available Slots:
#					- nsam:			number of samples for each iteration
#
#					- niter:		number of iteration
#
#					- theta:		mutation parameter
#
#					- obsVal:		vector with observed values for each test
#
#					- positions:	position of each polymorphic site
#
#					- trees:		if printtree = 1, find gene tree in Newick format, 
#									which is used in genetree() function
#
#					- seeds:		random numbers used to generate samples
#
#					- halplotypes:	haplotypes in each iteration
#
#					- stats:		variety of test stats compiled a matrix
#
#					- locProbLess:	Prob. that simulated val. <= to observed val. P(Sim <= Obs)
#
#					- locProbEqual:	Prob. that simulated val = to  observed val. P(Sim = Obs)
#
#					- locValidIter:	number of valid iteration for each test
#
#					- quantiles:	13 quantiles for each test
#
#
# AUTHOR:			Niroshan Nadarajah <niroshan.nadarajah@uni-duesseldorf.de>
#
# LAST MODIFIED:	10/11/03
#
###############################################################################

calc_probabilities <- function(cloci , niter, csStats, obsVal, lData,npop,testNames,numTests) {
	
	# lData beinhaltet zu cloci die statistiken fuer jede Population
	ntest <- dim(csStats[[1]])[2]
  
  # Init
  init <- as.matrix(vector("list",npop))
  csProbLess  <- init
  csProbEqual <- init
  csValidIter <- init
  quantiles   <- init
  ###########################
  popnames <- paste("pop",1:npop)
  colnames(csProbLess)  <- "cs.prob.less"
  colnames(csProbEqual) <- "cs.prob.equal"
  colnames(csValidIter) <- "cs.valid.iter"
  colnames(quantiles)   <- "quantiles"
  rownames(csProbLess)  <- popnames
  rownames(csProbEqual) <- popnames
  rownames(csValidIter) <- popnames
  rownames(quantiles)   <- popnames
  #print(csProbLess)

for(xx in 1:npop){	

  # get the number of tests 
	# ntest <- dim(csStats)[2]-1	
	# set up matrix and name columns for inserting the calculated quantiles
	quantiles[[xx]] <- matrix(nrow=length(quantileProbs), ncol=ntest)
	colnames(quantiles[[xx]]) <- testNames
	rownames(quantiles[[xx]]) <- quantileProbs
	
	# define space to save probability value for each locus
	csProbLess[[xx]]  <- matrix( ncol=ntest )
	csProbEqual[[xx]] <- matrix( ncol=ntest )
	csValidIter[[xx]] <- matrix( ncol=ntest )
	colnames(csProbLess[[xx]])  <- testNames
	colnames(csProbEqual[[xx]]) <- testNames
	colnames(csValidIter[[xx]]) <- testNames
	rownames(csProbLess[[xx]])  <- popnames[xx]
	rownames(csProbEqual[[xx]]) <- popnames[xx]
	rownames(csValidIter[[xx]]) <- popnames[xx]
	
	for (ctest in 1:ntest) {	
		
		# probability values can only be gained if observed values are present
		
    if (!is.na(obsVal[1][1]) ) {
		 	
			count 		  <- 0
			countequal 	  <- 0
			total 		  <- 0
			
			for (j in 1: niter) {
				currStat <- as.numeric(csStats[[xx]][,ctest][j])
		     
				if (!is.na(currStat) && !is.nan(currStat) && !is.na(obsVal[[xx]][,ctest])) {
					total <- total +1
					# probability that simulated values be smaller or equal than the observed value. p(sim <= obs): accuracy of +- 1e-6
					if ((obsVal[[xx]][,ctest] >= currStat+1e-05)  || (obsVal[[xx]][,ctest]  >= currStat -1e-05) ) {
						count = count + 1;
					}
					
					# probability that simulated values be equal to the observed value. p(sim = obs): accuracy of +- 1e-6
					if ( (obsVal[[xx]][,ctest] <= currStat+1e-05)  &&  (obsVal[[xx]][,ctest] >= currStat-1e-05) ) {
						countequal = countequal + 1;
					}
					
				}
			} #for loop j : niter
			
			#probability
			if ( total > 0 ) {
				csProbLess[[xx]]  [ctest]  <- count/total
				csProbEqual[[xx]] [ctest] <- countequal/total
				csValidIter[[xx]] [ctest] <- total
			}
		}
		
		# calculate values for specified quantiles
		qcol = calc_quantile(csStats[[xx]], ctest)
		quantiles[[xx]][,ctest] <- cbind(qcol)
		
	} #for loop ctest : ntest (number of tests performed)
 }# End for npops	
	# add all calculated values to the appropriate slots
	lData@quantiles   <- as.matrix(quantiles)
	
	
	if (length(obsVal) > 0) {
		
    #colnames(csProbLess) <- testNames
		lData@loc.prob.less  = as.matrix(csProbLess)
		
		#colnames(csProbEqual) <- testNames
		lData@loc.prob.equal = as.matrix(csProbEqual)
		
		#colnames(csValidIter) <- testNames
		lData@loc.valid.iter = as.matrix(csValidIter)
	}
	
	return(lData)
}

#------------------------------------------------------------------------------#
#								calc_quantile								   #
#------------------------------------------------------------------------------#
## this function returns the calculated quantiles for the column specified in 
## ctest. the quantile probabilties are specified in the vector quantileProbs
## which is found in the coalsim file
calc_quantile <- function(csStats, ctest) {
	return(quantile(csStats[,ctest], probs = quantileProbs, na.rm = TRUE))
}

