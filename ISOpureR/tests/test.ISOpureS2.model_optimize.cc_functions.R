# test.ISOpureS2.model_optimize.cc_functions.R ####################################################
# Testing script for functions needed in the optimilzation of mm in Step 2
# Test cc separately since it takes a while...
#
# The functions are of the form ISOpureS2.model_optimize.---
# (cc)
#   cc_loglikelihood
#   cc_deriv_loglikelihood

# LOAD DATA #######################################################################################
# load library 
library(ISOpureR);

# load the data from that path
data.path <-  paste0(file.path(system.file(package = "ISOpureR"), 'extdata', 'Beer'));  
load(file.path(data.path , 'beer.tumordata.1000.transcripts.30.patients.RData'));
load(file.path(data.path, 'beer.ISOpureS2model.1000.transcripts.30.patients.rounded.RData'));

# the normaldata and tumourdata should be matrices
beer.tumordata <- as.matrix(beer.tumordata);

# TEST CC FUNCTIONS ###############################################################################

# inputs needed for cc functions
# initial value of cc (this is testing just patient 1)
ww <- t(ISOpureS2model$cc_weights[1, , drop=FALSE]);

# test cc functions (just for patient 1, that's the third entry)
ISOpureS2.model_optimize.cc.cc_loglikelihood(ww, beer.tumordata, 1, ISOpureS2model);
ISOpureS2.model_optimize.cc.cc_deriv_loglikelihood(ww, beer.tumordata, 1, ISOpureS2model);

# remove parameters used in the test
rm(K, ww);