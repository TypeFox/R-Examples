### Generate phi.Obs based on Phi.

### log(Phi^{obs}) ~ N(log(Phi), sigma_W^2)
simu.phi.Obs <- function(Phi, sigmaW.lim = 1, bias.Phi = 0){
  tl.x <- length(Phi)
  orf.names <- names(Phi)

  ### Check orf.names.
  if(length(orf.names) != tl.x){
    orf.names <- paste("ORF", 1:tl.x, sep = "")
  }

  Phi <- log(Phi)
  Phi.lim <- range(Phi)

  if(length(sigmaW.lim) == 1){
    sigmaW.lim <- rep(sigmaW.lim, 2)
  }
  sigmaW.lim <- sigmaW.lim[1:2]
  sigmaW <- sigmaW.lim[2] -
            (Phi - Phi.lim[1]) / diff(Phi.lim) * diff(sigmaW.lim)

  error <- lapply(1:tl.x, function(i) rnorm(1, 0, sd = sigmaW[i]))
  error <- do.call("c", error)

  phi.Obs <- exp(Phi + error + bias.Phi)
  names(phi.Obs) <- orf.names

  phi.Obs
} # End of simu.phi.Obs().

