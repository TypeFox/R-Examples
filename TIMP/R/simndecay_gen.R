"simndecay_gen" <-
function (kinpar, tmax, deltat, specpar=vector(), lmin, lmax, deltal, 
          sigma, irf = FALSE, irfpar = vector(), seqmod = FALSE, dispmu
          = FALSE, nocolsums=FALSE, 
          disptau = FALSE, parmu = list(), partau = vector(), lambdac = 0, 
          fullk = FALSE, kmat = matrix(), jvec = vector(), specfun = "gaus", 
          nupow = 1, irffun = "gaus", kinscal = vector(), 
          lightregimespec = list(),
          specdisp = FALSE, specdisppar=list(), parmufunc = "exp",
          specdispindex = list(), amplitudes = vector(),specref=0, 
          fixedkmat=FALSE) 
{
    x <- seq(0, tmax, deltat)
    nt <- length(x)
    ncomp <- length(kinpar)
    x2 <- seq(lmin, lmax, deltal)
    nl <- length(x2)

    if(specdisp){
      ## store all the spectra; could do it otherwise if mem. is an issue
      EList <- list()
      for (i in 1:nt) {
        sp <- specparF(specpar = specpar, 
                            xi = x[i], 
                            i = i, specref = specref, 
                            specdispindex = specdispindex, 
                            specdisppar = specdisppar,
                            parmufunc = parmufunc)
        EList[[i]] <- calcEhiergaus(sp, x2, nupow)
      }
    }
    else 
      E2 <- calcEhiergaus(specpar, x2, nupow)
    
    if (!(dispmu || disptau)) {
      C2 <- compModel(k = kinpar, x = x, irfpar = irfpar, irf = irf, 
                      seqmod = seqmod, fullk = fullk, kmat = kmat,
                      jvec = jvec,amplitudes = amplitudes,
                      lightregimespec = lightregimespec,
                      nocolsums= nocolsums, kinscal = kinscal,
                      fixedkmat=fixedkmat)
      if(specdisp){
        psisim <- matrix(nrow = nt, ncol = nl)
        E2 <- EList[[1]]
        for (i in 1:nt) {
          psisim[i,] <- t(as.matrix(C2[i, ])) %*% t(EList[[i]])
        }
      }
      else 
        psisim <- C2 %*% t(E2)
    }
    else {
      psisim <- matrix(nrow = nt, ncol = nl)
      for (i in 1:nl) {
        irfvec <- irfparF(irfpar, lambdac, x2[i], i, dispmu, 
                          parmu, disptau, partau, "", "", "gaus")
        
        C2 <- compModel(k = kinpar, x = x, irfpar = irfpar, irf = irf, 
                        seqmod = seqmod, fullk = fullk, kmat = kmat,
                        jvec = jvec,amplitudes = amplitudes,
                        lightregimespec = lightregimespec,
                        nocolsums= nocolsums, kinscal = kinscal,
                      fixedkmat=fixedkmat)
        psisim[, i] <- C2 %*% cbind(E2[i, ])
      }
    }
    dim(psisim) <- c(nt * nl, 1)
    psi.df <- psisim + sigma * rnorm(nt * nl)
    dim(psi.df) <- c(nt, nl)
    
    kin(psi.df = psi.df, x = x, nt = nt, x2 = x2, nl = nl, C2 = C2, 
        E2 = E2, kinpar = kinpar, specpar = specpar, 
        seqmod = seqmod, irf = irf, irfpar = irfpar, 
        dispmu = dispmu, disptau = disptau, parmu = parmu, partau = partau, 
        lambdac = lambdac, simdata = TRUE, fullk = fullk, kmat = kmat, 
        jvec = jvec, fixedkmat=fixedkmat)
}

