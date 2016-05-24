# test.ISOpureS1.model_optimize.mm_functions.R #####################################
# Testing script for functions needed in the optimilzation of mm in Step 1
# Test mm separately since it takes a while...
#
# The functions are of the form ISOpureS1.model_optimize.---
# (mm)
#   mm_loglikelihood
#   mm_deriv_loglikelihood

# LOAD DATA #######################################################################################
# load library 
library(ISOpureR);

# load the data from that path
data.path <-  paste0(file.path(system.file(package = "ISOpureR"), 'extdata', 'Beer'));  
load(file.path(data.path , 'beer.tumordata.1000.transcripts.30.patients.RData'));
load(file.path(data.path, 'beer.ISOpureS1model.1000.transcripts.30.patients.RData'));

# the normaldata and tumourdata should be matrices
beer.tumordata <- as.matrix(beer.tumordata);

# TEST MM FUNCTIONS ###############################################################################

# inputs needed for mm functions
ww <- t(ISOpureS1model$mm_weights);

# test mm functions
ISOpureS1.model_optimize.mm.mm_loglikelihood(ww, beer.tumordata, ISOpureS1model);
ISOpureS1.model_optimize.mm.mm_deriv_loglikelihood(ww, beer.tumordata, ISOpureS1model);

# remove variables created for this test
rm(ww); 