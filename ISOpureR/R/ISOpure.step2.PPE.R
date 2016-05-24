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

### FUNCTION: ISOpure.step2.PPE.R ########################################################################
#
# This function performs the second step of the ISOpure purification algorithm, taking tumor data 
# normal profiles and returning the a list, ISOpureS2model, with all the updated parameters.
#
# Function call: ISOpure.step2.PPE <- function(tumordata, BB, ISOpureS1model, MIN_KAPPA)

### INPUT #########################################################################################
#  tumordata: (same as for ISOpureS1) a GxD matrix representing gene
#  expression profiles of heterogeneous (mixed) tumor samples, where G is
#  the number of genes, D is the number of tumor samples.
#
#  BB: (same as for ISOpureS1) represents B = [b_1 ... b_(K-1)] matrix (from
#      Genome Medicine paper). a Gx(K-1) matrix, where *K-1) is the number of normal profiles
#      (\beta_1,...,\beta_(K-1)), G is the number of genes.  These
#      are the normal profiles representing normal cells that contaminate
#      the tumor samples  (i.e. normal samples from the same tissue
#      location as the tumor).
#
#  ISOpureS1model: output model structure from ISOpureS1 code
#
#  MIN_KAPPA: (optional) The minimum value allowed for the strength parameters kappa_d placed
#  over the individual cancer profiles c_n (see Quon et al, 2013).  By default, this is set
#  to 1/min(m) (where m is the reference cancer profile) such that the log likelihood of the model 
#  is always finite.  However, when the min(m) is very small, this forces MIN_KAPPA to be very large, 
#  and can sometimes cause the cancer profiles to look too similar to the reference profile m
#  If this is the case, you can try setting MIN_KAPPA=1, or some other small value.  For reference, for the data
#  presented in Quon et al., 2013, MIN_KAPPA is on the order of 10^5 - 10^6.
#

### OUTPUT ########################################################################################
#
# ISOpureS2model: a list with the following important fields:
#
#  theta: a DxK matrix, giving the fractional composition of each tumor
#  sample.  Each row represents a tumor sample that was part of the input,
#  and the first K-1 columns correspond to the fractional composition with
#  respect to the Source Panel contaminants.  The last column represents
#  the fractional composition of the pure cancer cells.  In other words,
#  each row sums to 1, and element (i,j) of the matrix denotes the fraction 
#  of tumor i attributable to component j (where the last column refers to 
#  cancer cells, and the first K-1 columns refer to different 'normal cell'
#  components).  The "# cancer", or tumor purity, estimate of each tumor is
#  simply the last column of theta. 
#
#  alphapurities:  (same as ISOpureS1) tumor purities (alpha_i in paper),
#  same as the last column of the theta variable, pulled out for user 
#  convenience. 
#
#  cc_cancerprofiles: purified cancer profiles.  This matrix is of the same
#  dimensionality as tumordata, and is also on the same scale (i.e.
#  although ISOpureS2 treats purified cancer profiles as parameters of a
#  multinomial distribution, we re-scale them to be on the same scale as
#  the input tumor profiles -- see Genome Medicine paper).  column ii of
#  cc_cancerprofiles corresponds to column ii of tumordata.
# 
#  total_loglikelihood: log likelihood of the final model
# 
#  The model also includes several internal parameters: 
#    - the hyper-parameters vv and kappa from the Dirichlet distributions
#    - theta_weights and cc_weights, log_cc, which are used in the optimization
#      of theta and cc (instead of performing constrained optimization on these 
#      positively constrained variables directly, we optimize the log of them in an 
#      unconstrained fashion.) 
#    - log_BBtranspose are PPtranspose  are used in the calculations of loglikelihood 
#    - MIN_KAPPA as described above
#    - omega

ISOpure.step2.PPE <- function(tumordata, BB, ISOpureS1model, MIN_KAPPA=NULL, logging.level="INFO") {

	flog.threshold(logging.level);

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
		stop('Negative elements found in input matrix BB')
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

	# initial value of kappa of 10^4 seems to work well.  This is optimized later.
	kappa <- 10^4;

	# we work with transpose of BB
	BBtranspose <- t(BB);
	# initial starting point of theta will be where we left off at in ISOpureS1
	INITIAL_THETA <- ISOpureS1model$theta;
	
	# will start VV variable where we left off in ISOpureS1, except add some
	# small random amount to each component, because many components may be
	# exactly 1 (lower bound set on VV)
	INITIAL_VV <- ISOpureS1model$vv+runif(length(ISOpureS1model$vv));

	# prior on the tumor-specific cancer profiles is just the reference cancer
	# profile learned in ISOpureS1
	PPtranspose <- exp(ISOpureS1model$log_all_rates[nrow(ISOpureS1model$log_all_rates),]);

	NTOPICS <- dim(BBtranspose)[1]+1;

	# make sure PP and BB are scaled to sum to 1
	# rowsums changed to sum for PPtranspose since PPtranspose will be only one row 
	PPtranspose <-  PPtranspose / rep(sum(PPtranspose), nrow=1, ncol=ncol(PPtranspose));
	BBtranspose <- BBtranspose / ISOpure.util.repmat(rowSums(BBtranspose), 1, ncol(BBtranspose));

	total_loglikelihood <- -Inf;

	flog.info('Initializing model ...');

	# identify the minimum value of kappa such that the Dirichlet prior over
	# cancer profiles will give real-valued likelihoods

	if (is.null(MIN_KAPPA)){
		MIN_KAPPA <- 1/min(min(PPtranspose));
	}
	flog.info('MIN_KAPPA set to %s', MIN_KAPPA);

	# initialize the model structure that holds the parameters
	INIT_MODEL <- ISOpureS2.model_core.new_model(tumordata, max(kappa, 10*MIN_KAPPA), INITIAL_VV, PPtranspose, BBtranspose);
	INIT_MODEL$MIN_KAPPA <- MIN_KAPPA;
	INIT_MODEL$total_loglikelihood <- total_loglikelihood;

	# set initial theta to values learned in ISOpureS1
	INIT_MODEL$theta <- INITIAL_THETA;

	# estimate parameters/hidden variables
	ISOpureS2model <- ISOpureS2.model_core.optmodel(tumordata, INIT_MODEL);

	# copy over some important variables
	ISOpureS2model$alphapurities <- ISOpureS1model$alphapurities;
	
	# transpose log_cc to get same dimensions as tumordata, then scale
	# multinomial parameters to be on the same scale as tumordata
	ISOpureS2model$cc_cancerprofiles <- exp(t(ISOpureS2model$log_cc));
	ISOpureS2model$cc_cancerprofiles <- ISOpureS2model$cc_cancerprofiles * ISOpure.util.repmat(matrix(colSums(tumordata), nrow=1, ncol=ncol(tumordata)), nrow(ISOpureS2model$cc_cancerprofiles),1);
	
	return(ISOpureS2model);
}