MCMCMLE <- function(mc.num.iterations,
                    tolerance,
                    theta,
                    seed2 ,
					          possible.stats,
					          GERGM_Object,
					          force_x_theta_updates,
					          verbose) {

  # get MPLE thetas
  MPLE_Results <- run_mple(GERGM_Object = GERGM_Object,
                           verbose = verbose,
                           seed2 = seed2,
                           possible.stats = possible.stats)

  GERGM_Object <- MPLE_Results$GERGM_Object
  theta <- MPLE_Results$theta
  statistics <- MPLE_Results$statistics
  init.statistics <- MPLE_Results$init.statistics

  ##########################################################################
  ## Simulate new networks
  FIX_DEGENERACY <- FALSE
  for (i in 1:mc.num.iterations) {

    if (FIX_DEGENERACY) {
      MPLE_Results <- run_mple(GERGM_Object = GERGM_Object,
                               verbose = verbose,
                               seed2 = seed2,
                               possible.stats = possible.stats)

      GERGM_Object <- MPLE_Results$GERGM_Object
      theta <- MPLE_Results$theta
      statistics <- MPLE_Results$statistics
      init.statistics <- MPLE_Results$init.statistics
      FIX_DEGENERACY <- FALSE
    }
    GERGM_Object@theta.par <- as.numeric(theta$par)

    # now optimize the proposal variance if we are using Metropolis Hasings
    if (GERGM_Object@hyperparameter_optimization){
      if (GERGM_Object@estimation_method == "Metropolis") {
        GERGM_Object@proposal_variance <- Optimize_Proposal_Variance(
          GERGM_Object = GERGM_Object,
          seed2 = seed2,
          possible.stats = possible.stats,
          verbose = verbose)
        cat("Proposal variance optimization complete! Proposal variance is:",
            GERGM_Object@proposal_variance,"\n",
            "--------- END HYPERPARAMETER OPTIMIZATION ---------",
            "\n\n")
      }
    }

    GERGM_Object <- Simulate_GERGM(GERGM_Object,
                           seed1 = seed2,
                           possible.stats = possible.stats,
                           verbose = verbose)

    hsn <- GERGM_Object@MCMC_output$Statistics[,which(statistics == 1)]
    hsn.tot <- GERGM_Object@MCMC_output$Statistics

    # deal with case where we only have one statistic
    if(class(hsn.tot) == "numeric"){
      hsn.tot <- matrix(hsn.tot,ncol =1,nrow = length(hsn.tot))
      stats.data <- data.frame(Observed = init.statistics,
                               Simulated = mean(hsn.tot))
    }else{
      stats.data <- data.frame(Observed = init.statistics,
                               Simulated = colMeans(hsn.tot))
    }

    rownames(stats.data) <- possible.stats
    cat("Simulated (averages) and observed network statistics...\n")
    print(stats.data)
    GERGM_Object <- store_console_output(GERGM_Object,toString(stats.data))
    if(verbose){
      cat("\nOptimizing theta estimates... \n")
    }
    GERGM_Object <- store_console_output(GERGM_Object,"\nOptimizing Theta Estimates... \n")
    if(verbose){
    theta.new <- optim(par = theta$par,
                       log.l,
                       alpha = GERGM_Object@reduced_weights,
                       hsnet = hsn,
                       ltheta = as.numeric(theta$par),
                       together = GERGM_Object@downweight_statistics_together,
                       possible.stats= possible.stats,
                       GERGM_Object = GERGM_Object,
                       method = "BFGS",
                       hessian = T,
                       control = list(fnscale = -1, trace = 6))
    }else{
      theta.new <- optim(par = theta$par,
                         log.l,
                         alpha = GERGM_Object@reduced_weights,
                         hsnet = hsn,
                         ltheta = as.numeric(theta$par),
                         together = GERGM_Object@downweight_statistics_together,
                         possible.stats= possible.stats,
                         GERGM_Object = GERGM_Object,
                         method = "BFGS",
                         hessian = T,
                         control = list(fnscale = -1, trace = 0))
    }
    if(verbose){
      cat("\n", "Theta Estimates: ", paste0(theta.new$par,collapse = " "), "\n",sep = "")
    }
    GERGM_Object <- store_console_output(GERGM_Object,paste("\n", "Theta Estimates: ", paste0(theta.new$par,collapse = " "), "\n",sep = ""))
    theta.std.errors <- 1 / sqrt(abs(diag(theta.new$hessian)))
    # Calculate the p-value based on a z-test of differences
    # The tolerance is the alpha at which differences are significant
    p.value <- rep(0,length(as.numeric(theta$par)))
    count <- rep(0, length(as.numeric(theta$par)))
    for(j in 1:length(theta$par)){
      #two sided z test
      p.value[j] <- 2*pnorm(-abs((as.numeric(theta.new$par)[j] - as.numeric(theta$par)[j])/theta.std.errors[j]))
      #abs(theta.new$par[i] - theta$par[i]) > bounds[i]
      #if we reject any of the tests then convergence has not been reached!
      if(p.value[j] < tolerance){count[j] = 1}
    }
    if(verbose){
      cat("\np.values for two-sided z-test of difference between current and updated theta estimates:\n\n")
    }
    GERGM_Object <- store_console_output(GERGM_Object,"\np.values for two-sided z-test of difference between current and updated theta estimates:\n\n")
    if(verbose){
      cat(round(p.value,3), "\n \n")
    }
    GERGM_Object <- store_console_output(GERGM_Object,paste(p.value, "\n \n"))

    # calculate MCMC chain convergence diagnostic
    geweke_stat <- as.numeric(coda::geweke.diag(hsn.tot$edges)$z)
    cat("MCMC convergence Geweke test statistic:",geweke_stat,
        "\n(If the absolute value is greater than 1.7, increase MCMC_burnin)\n")
    GERGM_Object <- store_console_output(GERGM_Object,
      paste("MCMC convergence Geweke test statistic:",geweke_stat,
      "\n(If the absolute value is greater than 1.7, increase MCMC_burnin)\n"))

    allow_convergence <- TRUE
    if(max(abs(theta.new$par)) > 10000000){
      if(GERGM_Object@hyperparameter_optimization){
        message("Parameter estimates appear to have become degenerate, attempting to fix the problem...")
        GERGM_Object <- store_console_output(GERGM_Object,"Parameter estimates appear to have become degenerate, attempting to fix the problem...")
        # do not allow convergence
        allow_convergence <- FALSE
        # If we are using Metropolis Hastings, then try reducing weights and
        # upping the gain factor
        if (GERGM_Object@estimation_method == "Metropolis") {
          GERGM_Object@weights <- 0.9 * GERGM_Object@weights
          cat("Reducing exponential weights by 10 percent to:",
              GERGM_Object@weights,
              "in an attempt to address degeneracy issue...\n")
          GERGM_Object <- store_console_output(GERGM_Object,paste(
            "Reducing exponential weights by 10 percent to:",
            GERGM_Object@weights,
            "in an attempt to address degeneracy issue..."))

          if (GERGM_Object@MPLE_gain_factor == 0) {
            GERGM_Object@MPLE_gain_factor <- 0.05
          } else {
            GERGM_Object@MPLE_gain_factor <- 1.05 * GERGM_Object@MPLE_gain_factor
          }
          cat("Increasing MPLE gain factor by 5 percent percent to:",
              GERGM_Object@MPLE_gain_factor,
              "as exponential weights have decreased...\n")
          GERGM_Object <- store_console_output(GERGM_Object,paste(
            "Increasing MPLE gain factor by 5 percent percent to:",
            GERGM_Object@MPLE_gain_factor,
            "as exponential weights have decreased..."))
          # re-estimate thetas with more downweighting
          FIX_DEGENERACY <- TRUE
        }
        # additionally, try doubling the burin and the number of MCMC iterations.
        old_nsim <- GERGM_Object@number_of_simulations
        old_burinin <- GERGM_Object@burnin
        new_nsim <- 2 * old_nsim
        new_burnin <- 2 * old_burinin
        GERGM_Object@number_of_simulations <- new_nsim
        GERGM_Object@burnin <- new_burnin
        cat("Doubling burnin from:", old_burinin, "to", new_burnin,
            "and number of networks simulated from:", old_nsim, "to", new_nsim,
            "in an attempt to address degeneracy issue...\n")
        GERGM_Object <- store_console_output(GERGM_Object,paste(
          "Doubling burnin from:", old_burinin, "to", new_burnin,
          "and number of networks simulated from:", old_nsim, "to", new_nsim,
          "in an attempt to address degeneracy issue..."))
      }else{
        message("Parameter estimates appear to have become degenerate, returning previous thetas. Model output should not be trusted. Try specifying a larger number of simulations or a different parameterization.")
        GERGM_Object <- store_console_output(GERGM_Object,"Parameter estimates appear to have become degenerate, returning previous thetas. Model output should not be trusted. Try specifying a larger number of simulations or a different parameterization.")
        return(list(theta.new,GERGM_Object))
      }
    } else if (abs(geweke_stat) > 1.7){
      # if model was not degenerate but Geweke statistics say it did not converge
      # double number of iterations and burnin automatically.
      if (GERGM_Object@hyperparameter_optimization){
        old_nsim <- GERGM_Object@number_of_simulations
        old_burinin <- GERGM_Object@burnin
        new_nsim <- 2 * old_nsim
        new_burnin <- 2 * old_burinin
        GERGM_Object@number_of_simulations <- new_nsim
        GERGM_Object@burnin <- new_burnin
        cat("Doubling burnin from:", old_burinin, "to", new_burnin,
            "and number of networks simulated from:", old_nsim, "to", new_nsim,
            "in an attempt to address degeneracy issue...\n")
        GERGM_Object <- store_console_output(GERGM_Object,paste(
          "Doubling burnin from:", old_burinin, "to", new_burnin,
          "and number of networks simulated from:", old_nsim, "to", new_nsim,
          "in an attempt to address degeneracy issue..."))
        # do not allow convergence
        allow_convergence <- FALSE
      }
    }



    # check to see if we had a zero percent accept rate if using MH, and if so,
    # then adjust proposal variance and try again -- do not signal convergence.
    if (GERGM_Object@estimation_method == "Metropolis") {
      if (GERGM_Object@MCMC_output$Acceptance.rate == 0){
        old <- GERGM_Object@proposal_variance
        new <- old/2
        cat("Acceptance rate was zero, decreasing proposal variance from",old,
            "to",new,"and simulating a new set of networks...\n")
        GERGM_Object@proposal_variance <- new
        allow_convergence <- FALSE
      }
    }
    if (allow_convergence) {
      if (sum(count) == 0){
        #conditional to check and see if we are requiring more updates
        if(i >= force_x_theta_updates){
          if(verbose){
            message("Parameter estimates have converged")
          }
          GERGM_Object <- store_console_output(GERGM_Object,
                            "Parameter estimates have converged")
          GERGM_Object@theta_estimation_converged <- TRUE
          return(list(theta.new,GERGM_Object))
        }else{
          if(verbose){
            message(paste("Forcing",force_x_theta_updates,
                          "iterations of theta updates..."),sep = " ")
          }
          GERGM_Object <- store_console_output(GERGM_Object,paste("Forcing",
                            force_x_theta_updates,
                            "iterations of theta updates..."))
        }
      }
      # only updat parameter estimates if we had an acceptance rate greater than zero
      theta <- theta.new
      GERGM_Object@theta.par <- as.numeric(theta$par)
    }
  } #loop over MCMC outer iterations
  return(list(theta.new,GERGM_Object))
}
