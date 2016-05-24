setMethod("plotter", signature(model="spec"),
          function (model,multimodel, multitheta, plotoptions) {  
            if(plotoptions@residplot) {
              plotResids(multimodel, multitheta, plotoptions) 
            }
            if(plotoptions@residtraces)
              plotTracesResids(multimodel, multitheta, plotoptions)
            if (!plotoptions@notraces) 
              plotSpectraSuper(multimodel, multitheta, plotoptions)
            if(plotoptions@writefit || plotoptions@writefitivo)
              writeFit(multimodel, multitheta, plotoptions)
            plotoptions@addest <- c("specpar")
            plotoptions@paropt <- par(mgp = c(2, 1, 0),
                                      mar=c(3,3,3,2), oma = c(1,0,4,0), 
                                      mfrow=c(plotoptions@summaryplotrow,
                                        plotoptions@summaryplotcol))
            
            if(dev.cur() != 1)
              dev.new()
            par(mgp = c(2, 1, 0),
                mar=c(3,3,3,2), oma = c(1,0,4,0), 
                mfrow=c(plotoptions@summaryplotrow,
                  plotoptions@summaryplotcol))
            nt <- model@nt
            nl <- model@nl
            x <- model@x
            x2 <- model@x2
            m <- multimodel@modellist
            t <- multitheta 
            resultlist <- multimodel@fit@resultlist
            conmax <- list()       
            conlist <- getSpecList(multimodel, t)

	    if(plotoptions@adddataimage){
	            for(i in 1:length(m)) {
			datarev <- t(apply(m[[i]]@psi.df, 2, rev))
			xr <- rev(m[[i]]@x)
			x2lab <- xr[seq(1,length(xr), length=10)]
			x2at <- m[[i]]@x[seq(1,length(xr), length=10)]
			image(x=m[[i]]@x2, y=m[[i]]@x, z=datarev, ylab = plotoptions@xlab, yaxt="n",main = "Data", xlab = plotoptions@ylab, col=tim.colors())
			axis(2, at=x2at, labels=x2lab)
            	}
	    }

            for(i in 1:length(m)) {
              contoplot <- conlist[[i]]
              matplot(m[[i]]@x, contoplot,
                      type = "l", 
                      add=!(i==1), lty=i, ylab = "concentration", 
                      xlab=plotoptions@xlab,main = "Concentrations")
            }
                                        # SPECTRA
            specList <- vector("list", length(m))
            for(i in 1:length(m)) {
              if(m[[i]]@timedep)
                specpar <- specparF(t[[i]]@specpar, m[[i]]@x[1], 
                                    1, m[[i]]@specref, m[[i]]@specdispindex, 
                                    t[[i]]@specdisppar, parmufunc = m[[i]]@parmufunc)
              else 
                specpar <- t[[i]]@specpar 
              
              spectra <- doClpConstr(specModel(specpar, m[[i]]),
                                     clp_ind = 1, 
                                     clpCon = m[[i]]@clpCon, clpequ = t[[i]]@clpequ, 
                                     num_clpequ = length(m[[i]]@clpequspec), 
                                     usecompnames0 = m[[i]]@usecompnames0, 
                                     usecompnamesequ = m[[i]]@usecompnamesequ)
              specList[[i]] <- spectra 
              matplot(m[[i]]@x2, spectra, type = "l", 
		      main = "Spectra", 
		      xlab = plotoptions@ylab, ylab="amplitude", lty = i, 
		      add = !(i ==1)) 
            }
            for(i in 1:length(m)) {
              matplot(m[[i]]@x2, normdat(specList[[i]]), type = "l", 
		      main = "Normalized spectra", xlab = plotoptions@ylab, 
		      ylab="amplitude", lty = i, add = !(i ==1)) 
            }
                                        # RESIDUALS 
                                        # make a list of resid matrix
            residlist <- svdresidlist <- list() 
            for(i in 1:length(m)) {
		residuals <- matrix(nrow = m[[i]]@nt, 
				    ncol = m[[i]]@nl)
		for(j in 1:length(resultlist[[i]]@resid)){ 
                  residuals[j,] <- resultlist[[i]]@resid[[j]]
		}
		svdresidlist[[length(svdresidlist)+1]] <- doSVD(residuals,2,2) 
		residlist[[length(residlist)+1]] <- residuals 
	}
	limd<- max(  max(residlist[[1]]), abs(min(residlist[[1]]))) 

        if (! (any(diff(x) <= 0) || any(diff(x2) <= 0)))
                image(x, x2, residlist[[1]], xlab = plotoptions@xlab, 
		ylab = plotoptions@ylab, main = "Residuals Dataset 1", 
		zlim=c(-limd,limd), col=diverge_hcl(40, h = c(0, 120), c = 60, 
                                      l = c(45, 90), power = 1.2))
            if(nt > 1 && nl > 1){ 
	      matplot(x, svdresidlist[[1]]$left, type = "l",
                      main = "Left sing. vectors residuals ", log = "x", xlab =
                      plotoptions@xlab, col = 1, ylab=plotoptions@ylab)
	      if(length(m) > 1){
                for(i in 2:length(m)) {
                  matlines(m[[i]]@x, 
                           svdresidlist[[i]]$left, log = "x", 
                           type ="l", col = i)
                }
	      }
	      matplot(x2, t(svdresidlist[[1]]$right), type = "l",  
                      main = "Right sing. vectors residuals ", xlab = plotoptions@xlab, 
	      col = 1, ylab=plotoptions@ylab)
	      if(length(m) > 1){
                for(i in 2:length(m)) {
                  matlines(m[[i]]@x2, 
                           t(svdresidlist[[i]]$right), type ="l", col = i)
                 }
	      }
	      plot(1:length(svdresidlist[[1]]$values),
                   log10(svdresidlist[[1]]$values), ylab=plotoptions@ylab, 
                   main = "Sing. values residuals", type= "b", xlab="")
	      if(length(m) > 1){
                for(i in 2:length(m)) {
                  lines(1:length(svdresidlist[[i]]$values),
                        log10(svdresidlist[[i]]$values), 
                        type = "b", col = i)
                }
	      }
            }
                                        # DATA
                                        # make a list of svd data	
            
            svddatalist <- list() 
            for(i in 1:length(m)) {
              svddatalist[[length(svddatalist)+1]] <- doSVD(
                                                            multimodel@data[[i]]@psi.df,2,2) 
            }
            if(nt > 1 && nl > 1){ 
	      matplot(x, svddatalist[[1]]$left, type = "l",
                      main = "Left sing. vectors data", log = "x", xlab = plotoptions@xlab, 
                      col = 1, ylab=plotoptions@ylab)
	      if(length(m) > 1){
                for(i in 2:length(m)) {
                  matlines(m[[i]]@x, 
                           svddatalist[[i]]$left, log = "x", 
                           type ="l", col = i )
                }
	      }
	      matplot(x2, t(svddatalist[[1]]$right), type = "l", ylab="", 
                      main = "Right sing. vectors data", xlab = plotoptions@ylab, col = 1)
	      if(length(m) > 1){
                for(i in 2:length(m)) {
                  matlines(m[[i]]@x2, 
                           t(svddatalist[[i]]$right), type ="l", col = i)
                }
	      }
	      plot(1:length(svddatalist[[1]]$values),
                   log10(svddatalist[[1]]$values), ylab ="",
                   main = "Sing. values data", type= "b", xlab="")
	      if(length(m) > 1){
                for(i in 2:length(m)) {
                  lines(1:length(svddatalist[[i]]$values),
                        log10(svddatalist[[i]]$values), 
                        type = "b", col = i)
                }
	      }
            }
                                        # PLOT DSCALING PER CLP
            
            for(i in 1:length(m)) {
	      pl <- FALSE
    	      if(length(m[[i]]@dscalspec$perclp) != 0 )
                if(m[[i]]@dscalspec$perclp){
		  if(!pl) {
		    plot(m[[i]]@x2, t[[i]]@drel, 
                         main = "Dataset scaling per clp", xlab= plotoptions@ylab, 
                         ylab="", type = "l")  
		    pl <- TRUE 
                  } 
	          else 
		    lines(m[[i]]@x2, t[[i]]@drel, type = "l", col = i) 
                }	        
        }
            if(length(plotoptions@title) != 0){
              tit <- plotoptions@title
              if(plotoptions@addfilename) tit <- paste(tit,m[[i]]@datafile)
            }
    else {
      tit <- ""
      if(plotoptions@addfilename) tit <- paste(tit, m[[i]]@datafile)
    }
            mtext(tit, side = 3, outer = TRUE, line = 1)
            par(las = 2)
                                        # MAKE PS
            if(dev.interactive() && length(plotoptions@makeps) != 0) {
              if(plotoptions@output == "pdf")
                pdev <- pdf 
              else  pdev <- postscript
              dev.print(device=pdev, 
		file=paste(plotoptions@makeps, 
                  "_summary.ps", plotoptions@output, sep=""))
            }
            par(mfrow=c(plotoptions@summaryplotrow,1), new=TRUE)
            writeEst(multimodel, multitheta, plotoptions)
            # TODO: replace functionality provided by displayEst
            displayEst(plotoptions)
            
            if(plotoptions@plotkinspec) {
              plotClp(multimodel, t, plotoptions) 
            }
          })	

