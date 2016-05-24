simmr_elicit <-
function(x,proportion_means=rep(0,x$n_sources),proportion_sds=rep(0.1,x$n_sources),n_sims=1000) {
  # x must be a simmr_input xect
  if(class(x)!='simmr_input') stop("Input argument x must have come from simmr_load")
  # proportion_means must be a vector of length x$n_sources
  if(length(proportion_means)!=x$n_sources) stop('proportions_means must be of same length as the number of sources')
  # proportion_sds must be a vector of length x$n_sources
  if(length(proportion_sds)!=x$n_sources) stop('proportions_sds must be of same length as the number of sources')

  cat('Running elicitation optimisation routine...\n')

  # Perform an optimisation to match the standard deviations to a normal distribution
  clr_opt_mean = function(pars) {
    means = c(pars[(1:(x$n_sources-1))],-sum(pars[(1:(x$n_sources-1))]))
    # Generate n_sims observations from this normal distribution
    f = matrix(stats::rnorm(n_sims*x$n_sources,mean=means,sd=rep(1,x$n_sources)),ncol=x$n_sources,nrow=n_sims,byrow=TRUE)
    p = as.matrix(compositions::clrInv(f))
    return(sum((boot::logit(apply(p,2,'mean'))-boot::logit(proportion_means))^2))
  }

  opt_mean = stats::optim(c(rep(0,(x$n_sources-1))),clr_opt_mean,method="SANN")
  if(opt_mean$convergence!=0) {
    warning("Optimisation for means did not converge properly. Please either increase n_sims or adjust proportion_means and proportion_sds")
  } else {
    cat("Mean optimisation successful.\n")
  }
  best_mean = c(opt_mean$par,-sum(opt_mean$par))

  options(warn=-1)
  clr_opt_sd = function(pars) {
    sd = pars[1:(x$n_sources)]
    # Generate n_sims observations from this normal distribution
    f = matrix(stats::rnorm(n_sims*x$n_sources,mean=best_mean,sd=sd),ncol=x$n_sources,nrow=n_sims,byrow=TRUE)
    p = as.matrix(compositions::clrInv(f))
    return(sum((boot::logit(apply(p,2,'sd'))-boot::logit(proportion_sds))^2))
  }

  opt_sd = stats::optim(rep(1,(x$n_sources)),clr_opt_sd,method="SANN")
  if(opt_sd$convergence!=0) {
    warning("Optimisation for stand deviations did not converge properly. Please either increase n_sims or adjust proportion_means and proportion_sds")
  } else {
    cat("Standard deviation optimisation successful.\n")
  }
  best_sd = opt_sd$par
  options(warn=0)

  best_f = matrix(stats::rnorm(n_sims*x$n_sources,mean=best_mean,sd=best_sd),ncol=x$n_sources,nrow=n_sims,byrow=TRUE)
  best_p = as.matrix(compositions::clrInv(best_f))

  cat('Best fit estimates provide proportion means of:\n')
  cat(round(apply(best_p,2,'mean'),3))
  cat('\n')
  cat('... and best fit standard deviations of:\n')
  cat(round(apply(best_p,2,'sd'),3))
  cat('\n')

  return(list(mean=best_mean,sd=best_sd))

}
