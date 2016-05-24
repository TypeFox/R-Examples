"compModel" <-
  function (k, kinscal = vector(), x, irfpar = vector(), irf = FALSE,
            seqmod = FALSE, fullk = FALSE, kmat = matrix(), jvec =
              vector(), dscalspec = list(), drel = vector(), cohspec =
              list(), coh = vector(), oscspec=list(), oscpar = vector(), lamb = 1, dataset = 1, cohirf =
              vector(), mirf = FALSE, measured_irf = vector(), convalg =
              1, shiftmea = vector(), speckin2 = list(), usekin2 =
              FALSE, kinpar2 = vector(), kin2scal = vector(), reftau =
              0, anispec = list(calcani=FALSE), anipar = vector(),
            cohcol = vector(), amplitudes = vector(), streakT=0,
            streak=FALSE, doublegaus = FALSE, multiplegaus = FALSE, fixedkmat=FALSE, irffun
            = "gaus", kinscalspecial = list(), kinscalspecialspec =
              list(), lightregimespec = list(), nocolsums=FALSE,
            numericalintegration = FALSE, initialvals = vector(),
            reactantstoichiometrymatrix = vector(), 
            stoichiometrymatrix = vector(),lreturnA=FALSE)

{
    lightdiff <- length(lightregimespec) > 0
    if(lightdiff){
      c.temp <- getLightDiff(k, x, kinscal, kmat, jvec, fixedkmat,
                             kinscalspecial,
                             kinscalspecialspec, lightregimespec)
      return(c.temp)
    }
    if (fullk) {
      eig <- fullKF(k, kinscal, kmat, jvec, fixedkmat, kinscalspecial,
                    kinscalspecialspec, nocolsums)
      k <- -eig$values
      A <- eig$A
    }
    else if(seqmod)  
      A <- calcB(k) 
    ## modify k for anisotropy exp in time 
    if(anispec$calcani) {
      if (! (seqmod || fullk))
        A <- matrix(1, nrow = length(k), ncol = length(k)) 
      if(anispec$angle[dataset] != "MA") {
        k <- getAniK(k=k, dataset=dataset, ani=anispec, anipar=anipar)
        A <- getAniA(A=A, dataset=dataset, ani=anispec, anipar=anipar)
      }
    }
    if (irf) {
      if(length(shiftmea) != 0)
        shiftmea <- shiftmea[[1]]
      c.temp <- calcCirf(k=k, x=x, irfpar=irfpar, mirf=mirf, 
                         measured_irf=measured_irf, convalg=convalg,
                         shiftmea=shiftmea, 
                         lamb=lamb, reftau = reftau, doublegaus = doublegaus, 
                         streak=streak, streakT = streakT, irffun = irffun)
    }
    else if(numericalintegration)
      c.temp <- calcD(k, x, initialvals, 
                      reactantstoichiometrymatrix, stoichiometrymatrix)
    else c.temp <- calcC(k, x)
    
    ## now expand A to account for super ani
    if(anispec$calcani) {
      A <- getAniSuper(A = A, ani=anispec)
    }
    if(seqmod || fullk || anispec$calcani)
      c.temp <- c.temp %*% A 
    if(usekin2) {
      c.temp <- cbind( kin2scal[1] * c.temp, kin2scal[2] * 
                         compModel(k=kinpar2, fullk = speckin2$fullk,
                                   kmat = speckin2$kmat,
                                   kinscal = kinscal, jvec = speckin2$jvec,
                                   x=x, irfpar=irfpar, irf=irf,
                                   lamb=lamb, dataset=dataset, mirf = mirf,
                                   measured_irf = measured_irf,
                                   convalg = 1, shiftmea = shiftmea))
    }
    if(!is.null(cohspec$type) && cohspec$type != "") {
      c.temp <- cbind(c.temp, compCoh(irfpar, x, cohspec, coh, dataset,
                                      cohirf, mirf = mirf,
                                      measured_irf = measured_irf, convalg =
                                        convalg, shiftmea = shiftmea,
                                      lamb = lamb, ani=anispec,
                                      anipar=anipar, cohcol = cohcol))
    }
    if(!is.null(oscspec$type) && oscspec$type != "") {
      c.temp <- cbind(c.temp, compOsc(irfpar, x, oscspec, oscpar))
    }
    if (length(drel) != 0) {
      if (dscalspec$perclp) 
        c.temp <- drel[lamb] * c.temp
      else c.temp <- drel * c.temp
    }
    if (length(amplitudes) == ncol(c.temp)) {
      y <- matrix(0,nrow=length(amplitudes),ncol=length(amplitudes))
      diag(y) <- amplitudes
      c.temp <- c.temp %*% y
    }
    if (lreturnA)
      A
    else
      c.temp
  }

