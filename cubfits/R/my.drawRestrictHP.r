### Do the M-H on lognormal hyper-parameters.
my.drawRestrictHP <- function(proplist, list.Curr, phi.Curr){
  ### Compute probability ratio of likelihood.
  lpr <- sum(dlnorm(phi.Curr, meanlog = proplist$nu.Phi,
                    sdlog = proplist$sigma.Phi, log = TRUE)) -
         sum(dlnorm(phi.Curr, meanlog = list.Curr$nu.Phi,
                    sdlog = list.Curr$sigma.Phi, log = TRUE))

  ### The prior \pi(sigma.Phi^2) \propto 1.
  lpr <- lpr + 0
  ### Or, should the prior be \pi(sigma.Phi^2) \propto 1 / sigma.Phi^2 ?
  ### i.e. log \pi(sigma.Phi^2) \propto -2 * log(sigma.Phi)
  # lpr <- lpr - 2 * (proplist$sigma.Phi - list.Curr$sigma.Phi)

  ### log Acceptance probability.
  logAcceptProb <- lpr - proplist$lir

  ### Error handling -- interpreting NaN etc. as ~= 0.
  if(!is.finite(logAcceptProb)){
    warning("log acceptance probability not finite in hyperparam draw")
    logAcceptProb <- -Inf
  }

  ### Run MH acceptance rule.
  if(-rexp(1) < logAcceptProb){
    ret <- proplist
    ret$accept <- 1
  } else{
    ret <- list.Curr
    ret$accept <- 0
  }

  ### Return.
  ret
} # End of my.drawRestrictHP().
