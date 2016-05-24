make.ranges <-
function(emul, fix.betas) {

  cat('Obtaining emulator parameter ranges for optimization...\n')
  
# Preliminaries #!+
  yrange      <- max(emul$Y.mat) - min(emul$Y.mat) 
  par.min     <- apply(emul$Theta.mat, 2, min) #Min of each parameter 
  par.max     <- apply(emul$Theta.mat, 2, max) #Max of each parameter 
  parrange    <- par.max - par.min             #Ranges of each parameter
  X.mat.df    <- as.data.frame(emul$X.mat) 
  
# Rho #!+
  rho.lower   <- 0
  rho.upper   <- 0.9999999

# Kappa #!+
  kappa.lower <- (1E-5)*yrange
  kappa.upper <- 30*yrange

# Zeta #!+
  zeta.lower  <- (1E-5)*yrange 
  zeta.upper  <- 30*yrange

# Beta parameters (if needed) #w
# Estimated beta values are in beta.est$coefficients
  if (!fix.betas) { #If betas are estimated #w
     beta.est    <- lm(emul$Y.mat ~ 0 + ., data=X.mat.df, offset=NULL) 
     beta.upper  <- 10*abs(beta.est$coefficients) 
     beta.lower  <- -10*abs(beta.est$coefficients)
     # An exception if some betas are exactly 0 (almost impossible!)
     beta.0ind   <- beta.est$coefficients == 0 
     beta.lower[beta.0ind] <- -1E-6 
     beta.upper[beta.0ind] <- 1E-6
} else {#If betas are fixed!+
     beta.lower=NULL
     beta.upper=NULL
}

# Range parameters #!+
  phi.lower   <- parrange/1E4 
  phi.upper   <- parrange*10


# All parameters #!+
  all.lower   <- c(rho.lower, kappa.lower, zeta.lower, beta.lower, phi.lower)
  all.upper   <- c(rho.upper, kappa.upper, zeta.upper, beta.upper, phi.upper)

# Names for the ranges vectors #!+
  namevec   <- c("rho", "kappa", "zeta")
  if (!fix.betas) { 
      for (i in 1:length(beta.lower)) namevec <- c(namevec, "beta") 
  }
  for (i in 1:length(phi.lower)) namevec  <- c(namevec, "phi")
  names(all.lower) <- namevec 
  names(all.upper) <- namevec

# Output  #!+
  all.ranges  <- list(all.lower=all.lower, all.upper=all.upper)
  all.ranges
}
