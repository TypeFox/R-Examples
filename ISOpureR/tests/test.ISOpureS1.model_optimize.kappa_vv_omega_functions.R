# test.ISOpureS1.model_optimize.kappa_vv_omega_functions.R ########################################
# Testing script for functions needed in the optimilzation of kappa, vv, and omega in Step 1
# i.e. for the vector parameters and the scalar kappa
#
# I thought it would be better to test all of them at the same time, as the data needs to be loaded
# and that might be the most time-consuming step.
#
# The functions are of the form ISOpureS1.model_optimize.---
# (kappa)
#	kappa_loglikelihood.R,
#	kappa_deriv_loglikelihood.R,
#	kappa_compute_loglikelihood.R.
# (vv)
#   vv_compute_loglikelihood
#   ISOpure.model_optimize.vv.vv_loglikelihood
#	ISOpure.model_optimize.vv.vv_deriv_loglikelihood
# (omega)
#   omega_loglikelihood
#	omega_deriv_loglikelihood
#	omega_compute_loglikelihood

# LOAD DATA #######################################################################################
# load library 
library(ISOpureR);

# load the data from that path
data.path <-  paste0(file.path(system.file(package = "ISOpureR"), 'extdata', 'Beer'));  
load(file.path(data.path , 'beer.tumordata.1000.transcripts.30.patients.RData'));
load(file.path(data.path, 'beer.ISOpureS1model.1000.transcripts.30.patients.RData'));

# the normaldata and tumourdata should be matrices
beer.tumordata <- as.matrix(beer.tumordata);

# TEST KAPPA FUNCTIONS ############################################################################

# inputs needed for kappa functions
kappa <- as.numeric(ISOpureS1model$kappa); 
log_kappa <- log(ISOpureS1model$kappa - ISOpureS1model$MIN_KAPPA); 

# test kappa functions
ISOpureS1.model_optimize.kappa.kappa_loglikelihood(log_kappa, beer.tumordata, ISOpureS1model);
ISOpureS1.model_optimize.kappa.kappa_deriv_loglikelihood(log_kappa, beer.tumordata, ISOpureS1model);
ISOpureS1.model_optimize.kappa.kappa_compute_loglikelihood(kappa, beer.tumordata, ISOpureS1model);

# remove variables created for this test
rm(kappa, log_kappa);


# TEST VV FUNCTIONS ###############################################################################

# inputs needed for vv functions
# K = number of normal profiles+1
K <- length(ISOpureS1model$vv); 
# reshape vv into a Kx1 matrix
vv <- matrix(ISOpureS1model$vv,nrow=K, ncol=1);
# note that we don't directly optimize vvs because it has constraints
# (must be >=1 to guarantee real-valued likelihoods).
ww <- as.matrix(ISOpure.util.matlab_log(vv-1));
sum_log_theta <- matrix(colSums(ISOpure.util.matlab_log(ISOpureS1model$theta)), nrow=1, ncol=K);

# test vv functions (just checking patient 1 - last entry of input)
ISOpure.model_optimize.vv.vv_loglikelihood(ww, sum_log_theta, 1);
ISOpure.model_optimize.vv.vv_deriv_loglikelihood(ww, sum_log_theta, 1);
ISOpureS1.model_optimize.vv.vv_compute_loglikelihood(ISOpureS1model$vv, sum_log_theta, ncol(beer.tumordata));

# remove variables created for this test
rm(K, vv, ww, sum_log_theta);


# TEST OMEGA FUNCTIONS ############################################################################

# inputs needed for omega functions
ww <- ISOpureS1model$omega_weights;

# test kappa functions
ISOpureS1.model_optimize.omega.omega_loglikelihood(ww, beer.tumordata, ISOpureS1model);
ISOpureS1.model_optimize.omega.omega_deriv_loglikelihood(ww, beer.tumordata, ISOpureS1model);
ISOpureS1.model_optimize.omega.omega_compute_loglikelihood(ISOpureS1model$omega, beer.tumordata, ISOpureS1model);

# remove variables created for this test
rm(ww);
