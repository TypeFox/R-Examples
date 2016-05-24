#' @export
print.mlmc.test <- function(x, ...) {
  with(x, {
    cat("\n")
    cat("**********************************************************\n")
    cat("*** Convergence tests, kurtosis, telescoping sum check ***\n")
    cat("**********************************************************\n")
    cat("\n l   ave(Pf-Pc)    ave(Pf)   var(Pf-Pc)    var(Pf)")
    cat("    kurtosis     check \n-------------------------")
    cat("--------------------------------------------------\n")

    for(l in 0:L) {
      cat(sprintf("%2d   %8.4e  %8.4e  %8.4e  %8.4e  %8.4e  %8.4e \n",
                  l, del1[l+1], del2[l+1], var1[l+1], var2[l+1], kur1[l+1], chk1[l+1]))
    }

    if( kur1[length(kur1)] > 100.0 ) {
      cat(sprintf("\n WARNING: kurtosis on finest level = %f \n", kur1[length(kur1)]))
      cat(" indicates MLMC correction dominated by a few rare paths; \n")
      cat(" for (information on the connection to variance of sample variances,\n")
      cat(" see http://mathworld.wolfram.com/SampleVarianceDistribution.html\n\n")
    }

    if( max(chk1) > 1.0 ) {
      cat("\n WARNING: maximum consistency error = %f \n", max(chk1))
      cat(" indicates identity E[Pf-Pc] = E[Pf] - E[Pc] not satisfied \n\n")
    }

    cat("\n******************************************************\n")
    cat("*** Linear regression estimates of MLMC parameters ***\n")
    cat("******************************************************\n")
    cat(sprintf("\n alpha in %f  (exponent for (MLMC weak convergence)\n", alpha))
    cat(sprintf(" beta  in %f  (exponent for (MLMC variance) \n", beta))
    cat(sprintf(" gamma in %f  (exponent for (MLMC cost) \n", gamma))

    cat("\n")
    cat("***************************** \n")
    cat("*** MLMC complexity tests *** \n")
    cat("***************************** \n\n")
    cat("  eps       value   mlmc_cost   std_cost  savings     N_l \n")
    cat("--------------------------------------------------------- \n")

    for(i in 1:length(eps.v)) {
      cat(sprintf("%.4f  %.4e  %.3e  %.3e  %7.2f ",
                  eps.v[i], P[i], mlmc_cost[i], std_cost[i], std_cost[i]/mlmc_cost[i]))
      cat(sprintf("%9d", Nl[[i]]))
      cat("\n")
    }

    cat("\n")
  })
}
