# test.ISOpureS2.model_optimize.kappa_vv_functions.R ##############################################
# Testing script for functions needed in the optimilzation of kappa and vv in Step 2
# i.e. for the vector parameters (kappa is a vector in Step 2)
#
# The functions are of the form ISOpureS2.model_optimize.---
# (kappa)
#	kappa_loglikelihood.R,
#	kappa_deriv_loglikelihood.R,
#	kappa_compute_loglikelihood.R.
# (vv)
#   vv_compute_loglikelihood

# LOAD DATA #######################################################################################
# load library 
library(ISOpureR);

# load the data from that path
data.path <-  paste0(file.path(system.file(package = "ISOpureR"), 'extdata', 'Beer')); 
load(file.path(data.path , 'beer.tumordata.1000.transcripts.30.patients.RData'));
load(file.path(data.path, 'beer.ISOpureS2model.1000.transcripts.30.patients.rounded.RData'));

# the normaldata and tumourdata should be matrices
beer.tumordata <- as.matrix(beer.tumordata);

# TEST KAPPA FUNCTIONS ############################################################################

# inputs needed for kappa functions
kappa <- ISOpureS2model$kappa; 
log_kappa <- t(log(ISOpureS2model$kappa - ISOpureS2model$MIN_KAPPA)); 

# test kappa functions
ISOpureS2.model_optimize.kappa.kappa_loglikelihood(log_kappa, ISOpureS2model);
ISOpureS2.model_optimize.kappa.kappa_deriv_loglikelihood(log_kappa, ISOpureS2model);
ISOpureS2.model_optimize.kappa.kappa_compute_loglikelihood(kappa, ISOpureS2model);

# remove variables created for this test
rm(kappa, log_kappa);


# TEST VV FUNCTIONS ###############################################################################

# inputs needed for vv functions
# K = number of normal profiles+1
K <- length(ISOpureS2model$vv); 
# reshape vv into a Kx1 matrix
vv <- matrix(ISOpureS2model$vv,nrow=K, ncol=1);
# note that we don't directly optimize vvs because it has constraints
# (must be >=1 to guarantee real-valued likelihoods).
ww <- as.matrix(ISOpure.util.matlab_log(vv-1));
sum_log_theta <- matrix(colSums(ISOpure.util.matlab_log(ISOpureS2model$theta)), nrow=1, ncol=K);

# test vv functions
ISOpureS2.model_optimize.vv.vv_compute_loglikelihood(ISOpureS2model$vv, sum_log_theta, ncol(beer.tumordata));

# remove variables created for this test
rm(K, vv, ww, sum_log_theta);