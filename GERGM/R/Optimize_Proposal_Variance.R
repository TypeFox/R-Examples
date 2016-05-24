Optimize_Proposal_Variance <- function(GERGM_Object,
                                       seed2,
                                       possible.stats,
                                       verbose){

  # create a temporary GERGM object to use in proposal variance
  # optimization
  Opt_Prop_Var <- GERGM_Object
  Opt_Prop_Var@number_of_simulations <- max(GERGM_Object@number_of_simulations/10
                                            ,1000)
  Opt_Prop_Var@burnin <-  max(GERGM_Object@burnin/10
                              ,1000)
  if ((1/GERGM_Object@thin) >  (Opt_Prop_Var@number_of_simulations/10)) {
    Opt_Prop_Var@thin <- 10/Opt_Prop_Var@number_of_simulations
  }

  FOUND_ACCEPTABLE_PROP_VAR <- FALSE
  Acceptable_Proposal_Variance <- GERGM_Object@proposal_variance
  dampening_counter <- 1
  # find the optimal proposal variance
  while (!FOUND_ACCEPTABLE_PROP_VAR) {
    cat("--------- START HYPERPARAMETER OPTIMIZATION ---------",
        "\nSimulating",Opt_Prop_Var@number_of_simulations,
        "networks with proposal variance:", Opt_Prop_Var@proposal_variance,"\n")
    Opt_Prop_Var <- Simulate_GERGM(Opt_Prop_Var,
                                   seed1 = seed2,
                                   possible.stats = possible.stats,
                                   verbose = verbose)
    ar <- Opt_Prop_Var@MCMC_output$Acceptance.rate
    cat("Current acceptance rate:", ar,"\n")
    lb <- GERGM_Object@target_accept_rate - 0.05
    ub <- GERGM_Object@target_accept_rate + 0.05

    if (lb > ar) {
      change <- (1/(1 + dampening_counter)) * Opt_Prop_Var@proposal_variance
      Opt_Prop_Var@proposal_variance <- Opt_Prop_Var@proposal_variance - change
    } else if (ub < ar) {
      change <- (1/(1 + dampening_counter)) * (0.5 - Opt_Prop_Var@proposal_variance)
      Opt_Prop_Var@proposal_variance <- Opt_Prop_Var@proposal_variance + change
    } else {
      Acceptable_Proposal_Variance <- Opt_Prop_Var@proposal_variance
      FOUND_ACCEPTABLE_PROP_VAR <- TRUE
    }
    dampening_counter <-  dampening_counter + 1
    if (dampening_counter > 10) {
      cat("Stopping optimization, more iterations will likely not improve results...\n")
      FOUND_ACCEPTABLE_PROP_VAR <- TRUE
    }
  }

  return(Acceptable_Proposal_Variance)
}
