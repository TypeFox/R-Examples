emulator <-
function(mpars, moutput, par.reg, time.reg, kappa0, zeta0,
                  myrel.tol=NULL, twice=FALSE, fix.betas=TRUE)  {  

    # PRELIMINARIES #!+
    #===========================================================================
    # Check that par.reg has the correct number of elements #!+
    m.par      <- dim(mpars$par)[1] #!+
    if (length(par.reg) != m.par) {
       stop("***ERROR***: par.reg has wrong number of elements")
    }


    # INITIALIZE EMULATOR #!+
    #===============================================================================
    cat("Initializing the emulator...\n")
    # Initialize emulator. Initial emulator beta parameters ($beta.vec) are estimated
    # using multiple linear regression. 
    init.emul <- initialize.emul(mpars, moutput, par.reg, time.reg, kappa0, zeta0) #!+

    cat("\nInitial regression parameters:\n")
    cat(init.emul$beta.vec, '\n\n')

    # Evaluate initial emulator likelihood #!+
    init.pars <- make.parvec(init.emul, fix.betas) #!+
    n.par     <- length(init.emul$t.vec) 
    p.par     <- dim(init.emul$Theta.mat)[1]
    if (fix.betas) { #!+
      beta.vec <- init.emul$beta.vec
    } else {
      beta.vec <- NULL
    }

    init.lik  <- emul.lik(init.pars, init.emul$Y.mat, init.emul$X.mat, init.emul$t.vec,
                          init.emul$Theta.mat, n.par, p.par, fix.betas, NULL, NULL,
                          beta.vec) #!+
    cat("Initial emulator likelihood is: ", init.lik, "\n\n")

    # OPTIMIZE EMULATOR #!+
    #===============================================================================
    final.emul <- optimize.emul(init.emul, fix.betas, twice, myrel.tol ) #!+

    # OUTPUT #!+
    #=========
    final.emul
}
