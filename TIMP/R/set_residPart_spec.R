setMethod("residPart", signature(model="spec"), 
          function (model, group, multimodel, thetalist, clpindepX, finished, 
                    returnX, rawtheta) 
          {
            psi <- vector()
            spectra <- matrix()
            if(returnX) 
              thetalist <-  getThetaCl(rawtheta, multimodel)
            if(finished) 
              rlist <- list(irfvec=list())
            for(i in 1:length(group)) {
	      m <-  multimodel@modellist[[group[[i]][2]]]  
	      t <-  thetalist[[group[[i]][2]]] 
	      psi <- append(psi, m@psi.weight[group[[i]][1], ]) 
	      if(m@timedep){
                if(m@specdisp)
                  specpar <- specparF(specpar = t@specpar, 
                                      xi = m@x[group[[i]][1]], 
                                      i = group[[i]][1], specref = m@specref, 
                                      specdispindex = m@specdispindex, 
                                      specdisppar = t@specdisppar,
                                      parmufunc = m@parmufunc)
		else 
                  specpar <- t@specpar 
		spectra_i <- specModel(specpar, m)
		if(m@weight)
                  spectra_i <- weightNL(spectra_i,m,group[[i]][1])	
		if(ncol(spectra_i) != 0) {
                  spectra <- if(!identical(spectra, matrix()))
                    rbind(spectra, spectra_i) 
                  else 
                    spectra_i
	        }
              }
              else {
                if(identical(spectra, matrix()))
                  spectra <- clpindepX[[group[[i]][2]]]
                else 
                  spectra <- rbind(spectra, clpindepX[[group[[i]][2]]])
                
	      }
              
	      if(finished) {
		rlist$irfvec[[group[[i]][1]]] <- 
                  if(m@timedep) 
                    specpar
                  else c(0,0)
	      }
            }
            ## apply constraints to clp, using spec. from 1st dataset in group
            m <- multimodel@modellist[[group[[1]][2]]]
            if(m@timedep) 
              spectra <- doClpConstr(spectra, clp_ind = group[[1]][1], 
                                     clpCon = m@clpCon, clpequ = t@clpequ, 
                                     num_clpequ = length(m@clpequspec), 
                                     usecompnames0 = m@usecompnames0, 
                                     usecompnamesequ = m@usecompnamesequ)

            
            retval <- getResidRet(spectra, psi, rlist, returnX, finished,
                                  multimodel@nnls, multimodel@algorithm) 
            retval
            
          }) 
