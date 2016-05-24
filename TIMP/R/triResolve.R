"triResolve" <- function(currModel, currTheta){
  ##plot the alpha results
  opt <- currModel@optlist[[1]]
  if (opt@plot) 
        plotter(currModel@modellist[[1]], currModel, currTheta, opt)
  m <- currModel@modellist
  specList <- getSpecList(currModel,currTheta)
  aList <- lapply(specList, t)
  newdata <- aList[[1]]
  for(i in 2:length(aList)) 
    newdata <- rbind(newdata, aList[[i]])

  plotTri(newdata)

  startres <- getStartTri(newdata,ncomp=nrow(aList[[1]]),nexp=length(aList))

  ampList <- startres$ampList
  fixedamps <- startres$fixedamps
   
  aDatList  <- list(dat(psi.df = newdata,
                        x = 1:nrow(newdata), nt = nrow(newdata),
                        x2 = 1:ncol(newdata), nl = ncol(newdata)))  
  
  amp_model <- initModel(mod_type = "amp", amps = ampList,
                         fixed = list(amps = fixedamps))
  
  fit <- fitModel(data = aDatList, modspec = list(amp_model),
                  opt = opt(iter=3, nnls = 2, plot=FALSE))
 
  newcp <- fit$currModel@fit@resultlist[[1]]@cp
  ## add the new amplitude estimates 
  for(i in 1:length(currModel@data)) {
    currTheta[[i]]@amplitudes <- fit$currTheta[[1]]@amps[[i]]
    currModel@fit@resultlist[[1]]@cp[[i]] <- newcp[[i]]
  }
  ## now the residuals have changed, and so has the fit
  ## need to calculate the new model
  ## but will just plug in the new estimates into a call to 
  ## fitModel (0 iter.) separately
	
  list(currModel=fit$currModel,currTheta=fit$currTheta) 
  
}



