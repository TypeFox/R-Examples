Convert_Simulated_Networks_To_Observed_Scale <- function(
  GERGM_Object){
  # determine the number of MCMC samples
  samples <- dim(GERGM_Object@MCMC_output$Networks)[3]
  num.nodes <- GERGM_Object@num_nodes
  triples = t(combn(1:num.nodes, 3))
  stats <- rep(1,length(GERGM_Object@stats_to_use))
  transformation_type <- GERGM_Object@transformation_type

  # figure out what kind of transformation we did, then convert back.
  printseq <- round(seq(1,samples, length.out = 11)[2:11],0)
  printcounter <- 1
  if (length(GERGM_Object@data_transformation) > 0) {
    if(GERGM_Object@is_correlation_network){
      cat("Currently not implemented for correlation networks with covariates.")
    }else{
      for(i in 1:samples){
        if(i == printseq[printcounter]){
          cat(10*printcounter,"% complete...\n", sep = "")
          printcounter <- printcounter +1
        }
        # if we did a transformation (which is the default if we are including an intercept)
        if(transformation_type == "logcauchy" | transformation_type == "cauchy"){
          GERGM_Object@MCMC_output$Networks[,,i] <- qst(
            GERGM_Object@MCMC_output$Networks[,,i],
            GERGM_Object@BZ,
            GERGM_Object@BZstdev,
            1)
          if(transformation_type == "logcauchy"){
            GERGM_Object@MCMC_output$Networks[,,i] <- exp(GERGM_Object@MCMC_output$Networks[,,i])
          }
          diag(GERGM_Object@MCMC_output$Networks[,,i]) <- 0
        }
        if(transformation_type == "lognormal" | transformation_type == "gaussian"){
          GERGM_Object@MCMC_output$Networks[,,i] <- qst(
            GERGM_Object@MCMC_output$Networks[,,i],
            GERGM_Object@BZ,
            GERGM_Object@BZstdev,
            Inf)
          if(transformation_type == "lognormal"){
            GERGM_Object@MCMC_output$Networks[,,i] <- exp(GERGM_Object@MCMC_output$Networks[,,i])
          }
          diag(GERGM_Object@MCMC_output$Networks[,,i]) <- 0
        }
        GERGM_Object@MCMC_output$Statistics[i,] <- h2(
          net = GERGM_Object@MCMC_output$Networks[,,i],
          triples = triples,
          statistics = stats,
          alphas = NULL,
          together = 1,
          directed = TRUE)
      }
    }

  }else{
    if(GERGM_Object@is_correlation_network){
      for(i in 1:samples){
        if(i == printseq[printcounter]){
          cat(10*printcounter,"% complete...\n", sep = "")
          printcounter <- printcounter +1
        }
        # symmetrize incase there were any numerical imperfections
        symnet <- Symmetrize_Network(GERGM_Object@MCMC_output$Networks[,,i])
        GERGM_Object@MCMC_output$Networks[,,i] <- bounded.to.correlations(
          symnet)
      }
    }else{
      # if we did not do a transformation (only structural terms)
      cat("Currently not implemented for non-transformed networks.")
    }
  }
  return(GERGM_Object)
}
