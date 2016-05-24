    setMethod("getClpindepX", signature(model = "kin"), function(model, 
        multimodel, theta, returnX, rawtheta, dind) {
	if(returnX) 
		 theta <-  getThetaCl(rawtheta, multimodel)[[dind]]
        x <- compModel(k = theta@kinpar, kinscal = theta@kinscal, 
            x = model@x, irfpar = theta@irfpar, irf = model@irf, 
	    seqmod = model@seqmod, 
            convalg = model@convalg, fullk = model@fullk, kmat = model@kmat, 
            jvec = model@jvec, dscalspec = model@dscalspec, drel = theta@drel, 
            mirf = model@mirf, measured_irf = model@measured_irf, 
	    speckin2 = model@speckin2, 
	    usekin2 = model@usekin2, kinpar2 = theta@kinpar2, 
	    kin2scal = theta@kin2scal, reftau = model@reftau, 
	    anispec = model@anispec, anipar = theta@anipar, 
	    cohcol = model@cohcol, amplitudes = theta@amplitudes, 
            streak = model@streak, streakT = model@streakT, 
	    doublegaus = model@doublegaus,fixedkmat=model@fixedkmat,
            irffun=model@irffun, kinscalspecial = theta@kinscalspecial,
            kinscalspecialspec = model@kinscalspecialspec,
                       lightregimespec = model@lightregimespec,
                       numericalintegration = model@numericalintegration,
          initialvals = model@initialvals,
          reactantstoichiometrymatrix = model@reactantstoichiometrymatrix,
          stoichiometrymatrix = model@stoichiometrymatrix)
	if(returnX) 
		    x <- as.vector(x) 
	
	x
		    
    })
