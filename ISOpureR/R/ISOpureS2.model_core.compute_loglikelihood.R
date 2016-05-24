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

### FUNCTION: ISOpureS2.model_core.compute_loglikelihood.R #############################################################
# 
# Input variables: 
#   tumordata: a GxD matrix representing gene expression profiles of tumor samples
#   model: list containing all the parameters updated in ISOpure step two iterations
#
# Output variables: 
#   loglikelihood: the scalar value of the full loglikelihood function that is being optimized in 
#       ISOpure step two

ISOpureS2.model_core.compute_loglikelihood <- function(tumordata, model) { 

    loglikelihood <- 0;

    # D is the number of tumor samples
    D <- dim(tumordata)[2];
    # G is the number of transcripts/genes
    G <- dim(model$log_BBtranspose)[2];

    # add the loglikelihood of the thetas, log P(\theta_d|\vv)
    loglikelihood <- loglikelihood + D*(lgamma(sum(model$vv)) - sum(lgamma(model$vv))) + sum((model$vv-1) %*% t(ISOpure.util.matlab_log(model$theta)));

    kappaomegaPP <- model$omega %*% model$PPtranspose * ISOpure.util.repmat(t(model$kappa), 1, G);
 
    # add the loglikelihood of the observed tumor profiles t_i
    for (dd in 1:D) {
        log_all_rates <- rbind(model$log_BBtranspose, model$log_cc[dd,,drop=FALSE]);

        # add ln P(t_d|\theta,BB,cc_i)
        log_P_t_given_theta <- ISOpure.util.logsum(t(ISOpure.util.repmat(ISOpure.util.matlab_log(model$theta[dd,,drop=FALSE]), G, 1)) + log_all_rates, 1);

        loglikelihood <- loglikelihood + (log_P_t_given_theta %*% tumordata[,dd]);

        # add ln P(cc_i | \kappa, \omega, PP)
        loglikelihood <- loglikelihood + lgamma(sum(kappaomegaPP[dd,])) - sum(lgamma(kappaomegaPP[dd,])) + ((kappaomegaPP[dd,,drop=FALSE]-1) %*% t(model$log_cc[dd,,drop=FALSE]));
    }
    return(loglikelihood);    
}
