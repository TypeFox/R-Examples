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

### FUNCTION: ISOpureS2.model_optimize.opt_cc.R ############################################################################
# 
# Input variables: 
#   tumordata: a GxD matrix representing gene expression profiles of tumour samples
#   model: list containing all the parameters to be optimized
#   NUM_ITERATIONS_RMINIMIZE: minimum number of iteration that the minimization algorithm runs
#   iter: the iteration number
#   NUM_GRID_SEARCH_ITERATIONS: number of times to try restarting with different initial values
#
# Output variables: 
#   model: the model with cc_weights and log_cc updated
#
# Notes: The goal of this function is to optimize the tumor-specific cancer profiles.
#   Because cc is constrained (each cc_i are parameters of multinomial/discrete distribution), we 
#   don't directly optimize the likelihood function w.r.t. cc, but we perform change of variables 
#   to do unconstrained optimization.  We therefore store these unconstrained variables in the
#   field "cc_weights", and update these variables.

ISOpureS2.model_optimize.opt_cc <- function(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS){
    
    # D = number of patients/tumour samples
    D <- nrow(model$log_cc);
    # G = number of genes
    G <- ncol(model$log_BBtranspose); 

    # if cc_weights is not a field (i.e. for the first iteration), initialize
    # cc_weights to the last row of log_all_rates, which in the ISOpure case is just a linear
    # combination of all the normal profiles with equal weights
    if (!any(names(model) == "cc_weights")) {
        model$cc_weights <- model$log_cc; 
    }

    # note: cc_weights is updated separately for each patient, so may be able to parallelize this?
    for (dd in 1:D){
        init_xx <- t(model$cc_weights[dd, , drop=FALSE]);
        
        # perform the optimization
        retval <- ISOpure.model_optimize.cg_code.rminimize(init_xx, ISOpureS2.model_optimize.cc.cc_loglikelihood, ISOpureS2.model_optimize.cc.cc_deriv_loglikelihood, NUM_ITERATIONS_RMINIMIZE, tumordata=tumordata, dd, model=model);

        # convert optimized values back to mm 
        xx <- as.vector(retval[[1]]); 
        model$cc_weights[dd,] <- xx;

        # subtracting ISOpure.util.logsum(mm_weights) from mm_weights corresponds to dividing by total to make sure sum is 1 (once take exponent)
        model$log_cc[dd,] <- xx-ISOpure.util.logsum(as.matrix(xx, nrow=1, ncol=length(xx)),1);
    } 
    
    return(model);
}
