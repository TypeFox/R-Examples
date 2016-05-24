### Propose bias.Phi in log scale via random walk.
my.proposebiasPhi.RW_Norm <- function(bias.Phi.Curr,
    bias.Phi.DrawScale = .CF.CONF$bias.Phi.DrawScale,
    bias.Phi.DrawScale.prev = .CF.CONF$bias.Phi.DrawScale){
  ### Draw from proposcal.
  bias.Phi.New <- rnorm(1, mean = bias.Phi.Curr, sd = bias.Phi.DrawScale)

  ### Compute log ratio of prior.
  lir <- 0    # no Jacobin for no transformation
  if(bias.Phi.DrawScale != bias.Phi.DrawScale.prev){
    lir <- dnorm(bias.Phi.New, bias.Phi.Curr,
                 sd = bias.Phi.DrawScale, log = TRUE) -
           dnorm(bias.Phi.Curr, bias.Phi.New,
                 sd = bias.Phi.DrawScale.prev, log = TRUE)
  }

  ### Return.
  ret <- list(bias.Phi = as.numeric(bias.Phi.New),
              lir = lir)
  ret
} # End of my.proposebiasPhi.RW_Norm().
