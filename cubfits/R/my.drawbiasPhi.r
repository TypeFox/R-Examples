### Do the M-H on Bias Phi.
my.drawbiasPhi <- function(proplist, list.Curr, log.phi.Obs, log.phi.Curr,
    sigmaW.Curr){
  ### Compute probability ratio.
  ### Since this fact, we can use dnorm in this function
  ### x <- 1.5; m <- 2; s <- 3
  ### dnorm(log(phi), m, s, log = TRUE) - log(phi) ==
  ###   dlnorm(phi, m, s, log = TRUE)
  lpr <- sum(dnorm(log.phi.Obs, log.phi.Curr + proplist$bias.Phi,
                   sd = sigmaW.Curr, log = TRUE)) -
         sum(dnorm(log.phi.Obs, log.phi.Curr + list.Curr$bias.Phi,
                   sd = sigmaW.Curr, log = TRUE))

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
} # End of my.drawbiasPhi().
