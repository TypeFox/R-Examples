# The ISOpureR package is copyright (c) 2014 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

### FUNCTION: ISOpure.step1.CPE.R ########################################################################
#
# This function performs the first step of the ISOpure purification algorithm, taking tumor data 
# normal profiles and returning the a list, ISOpureS1model, with all the updated parameters.
# 
# Function call: ISOpure.step1.CPE <- function(tumordata, BB, PP, MIN_KAPPA) 

### INPUT #########################################################################################
#  tumordata: a GxD matrix representing gene expression profiles of
#  heterogeneous (mixed) tumor samples, where G is the number of genes, D
#  is the number of tumor samples.
#
#  BB: represents B = [b_1 ... b_(K-1)] matrix (from Genome Medicine paper)
#      a Gx(K-1) matrix, where (K-1) is the number of normal profiles
#      (\beta_1,...,\beta_(K-1)), G is the number of genes.  These
#      are the normal profiles representing normal cells that contaminate
#      the tumor samples  (i.e. normal samples from the same tissue
#      location as the tumor).  The minimum element of BB must be greater than 0 --
#      i.e. every gene/transcript must be observed on some level in each normal sample.
#
#  PP: (optional) a GxM matrix, representing the expression profiles whose
#  convex combination form the prior over the purified cancer profile
#  learned.   If only primary tumors from the same site of origin are
#  represented in tumordata, then this will be ste to be the same as BB (default
#  behavior).   This parameter is for backwards compatibility and replacing
#  the original ISOLATE code from the 2009 Bioinformatics paper, and can
#  represent potential sites of origins of the metastatic tumor (in which
#  case tumordata represents one or more expression profiles of the
#  secondary tumor).  Set PP=BB for default behavior, and if you need to specify
#  MIN_KAPPA.
#
#  MIN_KAPPA: (optional) The minimum value allowed for the strength parameter kappa' placed
#  over the reference cancer profile m (see Quon et al, 2013).  By default, this is set
#  to 1/min(BB), such that the log likelihood of the model is always finite.  However,
#  when the min(BB) is very small, this forces MIN_KAPPA to be very large, and can sometimes
#  cause the reference profile m to look too much like a 'normal profile' (and therefore
#  you may observe the tumor samples having low % cancer content estimates).  If this is the case,
#  you can try setting MIN_KAPPA=1, or some other small value.  For reference, for the data
#  presented in Quon et al., 2013, MIN_KAPPA is on the order of 10^5. 

### OUTPUT ########################################################################################
#
# ISOpureS1model: a list with the following important fields:
#
#  theta: a DxK matrix, giving the fractional composition of each tumor
#  sample.  Each row represents a tumor sample that was part of the input,
#  and the first K-1 columns correspond to the fractional composition with
#  respect to the Source Panel contaminants.  The last column represents
#  the fractional composition of the pure cancer cells.  In other words,
#  each row sums to 1, and element (i,j) of the matrix denotes the fraction
#  of tumor i attributable to component j (where the last column refers to
#  cancer cells, and the first K-1 columns refer to different 'normal cell'
#  components).  The 'cancer', or tumor purity, estimate of each tumor is
#  simply the last column of theta.
#
#  alphapurities:  tumor purities (alpha_i in paper), same as the last
#  column of the theta variable, pulled out for user convenience.
#
#  mm: reference cancer profile, in the form of parameters of a multinomial
#  or discrete distribution (sum of elements is 1).  This is the same as
#  the purified cancer profile that ISOLATE was designed to learn.
#
#  omega: a Mx1 vector describing the convex combination weights learned by
#  ISOpure step 1 over the PPtranspose matrix, that when applied to the
#  Site of Origin Panel, forms the prior over the reference cancer profile.
#  When ISOpure step 1 is used in a similar fashion to the ISOLATE
#  algorithm, entry ii indicates the "probability" that the normal profile
#  in the (ii)th column of PP is the site of origin of the secondary tumors
#  stored in tumordata.  
#
#  total_loglikelihood: log likelihood of the model
#
#  The model also includes several internal parameters: 
#    - the hyper-parameters vv and kappa from the Dirichlet distributions
#    - mm_weights, theta_weights, and omega_weights which are used in the optimization
#      of mm, theta, and omega (instead of performing constrained optimization on these 
#      positively constrained variables directly, we optimize the log of them in an 
#      unconstrained fashion.) 
#    - log_BBtranspose, PPtranspose, and log_all_rates are used in the calculations of 
#      loglikelihood 
#    - MIN_KAPPA as described above

ISOpure.step1.CPE <- function(tumordata, BB, PP=NULL, MIN_KAPPA=NULL, logging.level="INFO") { 

	flog.threshold(logging.level);

	# by default, we are looking at primary tumors from the same site, so the
	# "Source Panel" (profiles that form the components of the prior over the reference 
	# cancer profile) is the same as BB
	
	# I tried to set PP=BB in the function definition, but it didn't seem to work
	# thought this is a better check than using nargs() < 3 as for the Matlab as MIN_KAPPA
	# may be defined, but PP not defined 
	if (is.null(PP)) {
		PP <- BB;
	}

	# Tumordata --------------------------------------------------------------#

	# make sure tumordata is a proper intensity/read count matrix (no negative elements)
	if (min(min(tumordata))<0) {
		flog.fatal('Negative elements found in input matrix tumordata');
		stop('Negative elements found in input matrix tumordata');
	}

	# make sure minimum value is not 0 (i.e. all genes need to have some probability of being observed in a sample)
	if (min(min(tumordata))==0) { 
		nzix <- which(tumordata>0);
		mymin <- min(tumordata[nzix]);
		tumordata[which(tumordata==0)] <- mymin;
		flog.warn('Minimum element in input matrix tumordata is 0 -- setting all zeros to smallest non-zero element, %s', mymin); 
	}

	# make sure data is not log transformed
	if (max(max(tumordata)) < 30) {
		flog.warn('Maximum element in matrix tumordata is less than 30 -- make sure data is in normal (not log) space, otherwise output is wrong.');
	}

	# BB (Normaldata) --------------------------------------------------------#

	# make sure BB is a proper intensity/read count matrix (no negative elements)
	if (min(min(BB))<0) {
		flog.fatal('Negative elements found in input matrix BB');
		stop('Negative elements found in input matrix BB');
	}

	# make sure minimum value is not 0 (i.e. all genes need to have some probability of being observed in a sample)
	if (min(min(BB))==0) { 
		nzix <- which(BB>0);
		mymin <- min(BB[nzix]);
		BB[which(BB==0)] <- mymin;
		flog.warn('Minimum element in input matrix BB is 0 -- setting all zeros to smallest non-zero element, %s', mymin); 
	}

	# make sure data is not log transformed
	if (max(max(BB)) < 30) {
		flog.warn('Maximum element in matrix BB is less than 30 -- make sure data is in normal (not log) space, otherwise output is wrong.');
	}

	# PP (Normaldata) --------------------------------------------------------#

	# do same checks for the PP matrix
	if (min(min(PP))<0) {
		flog.fatal('Negative elements found in input matrix PP');
		stop('Negative elements found in input matrix PP');
	}

	# make sure minimum value is not 0 (i.e. all genes need to have some probability of being observed in a sample)
	if (min(min(PP))==0) { 
		nzix<-which(PP>0);
		mymin<-min(PP[nzix]);
		PP[which(PP==0)]<-mymin;
		flog.warn('Minimum element in input matrix PP is 0 -- setting all zeros to smallest non-zero element, %s', mymin);
	}

	# make sure data is not log transformed
	if (max(max(PP)) < 30) {
		flog.warn('Maximum element in matrix PP is less than 30 -- make sure data is in normal (not log) space, otherwise output is wrong.');
	}

	# we work with these BB and PP in their transpose
	BBtranspose <- t(BB);
	PPtranspose <- t(PP);

	# NTOPICS is the total number of component profiles (# normals + 1 for the
	# reference cancer profile)
	NTOPICS <- dim(BBtranspose)[1]+1;

	#initial value of kappa of 10^4 seems to work well.  This is optimized later.
	kappa <- 10^4;

	# also randomly initialize prior over mixing proportions, theta 
	INITIAL_VV <- runif(NTOPICS) + 1;     
	# last entry weighed more heavily
	INITIAL_VV[length(INITIAL_VV)] <- INITIAL_VV[length(INITIAL_VV)] + 5;
	
	# we standardize the normal profiles to sum to 1, so that we can interpret
	# them as parameters of a discrete or multinomial distribution
	PPtranspose <-  PPtranspose / ISOpure.util.repmat(rowSums(PPtranspose), 1, ncol(PPtranspose));
	BBtranspose <- BBtranspose / ISOpure.util.repmat(rowSums(BBtranspose), 1, ncol(BBtranspose));

	total_loglikelihood <- -Inf;

	flog.info('Initializing model...');

	# MIN_KAPPA represents the minimum value of kappa that we enforce during
	# optimization, so that the Dirichlet distribution gives real valued
	# loglikelihood values

	# again use the function is.null instead of nargs() to set default
	if (is.null(MIN_KAPPA)){
		MIN_KAPPA <- 1/min(min(PPtranspose));
	}
	flog.info('MIN_KAPPA set to %s', MIN_KAPPA);

	# initialize the model structure that holds the parameters 
	INIT_MODEL <- ISOpureS1.model_core.new_model(tumordata, max(kappa, 10*MIN_KAPPA), INITIAL_VV, PPtranspose, BBtranspose);
	INIT_MODEL$MIN_KAPPA <- MIN_KAPPA;
	INIT_MODEL$total_loglikelihood <- total_loglikelihood;

	ISOpureS1model <- ISOpureS1.model_core.optmodel(tumordata, INIT_MODEL);
	
	# explicitly save reference cancer profile mm, cancer purities alpha for
	# ease of extraction by user
	ISOpureS1model$mm <- exp(ISOpureS1model$log_all_rates[nrow(ISOpureS1model$log_all_rates),]);
	ISOpureS1model$alphapurities <- ISOpureS1model$theta[,ncol(ISOpureS1model$theta)];

	return(ISOpureS1model);
}