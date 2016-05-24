gridfn <- function(varname, gridsteps=blackbox.getOption("gridstepsNbr"), margefrac=0, varnameS=NULL) {
  fittedNames <- blackbox.getOption("fittedNames")
  plotRange <- blackbox.getOption("plotRange")
  Grid <- c();
  if ( ! is.null(plotRange[[varname]]) ) { ## plotRange as given by user, overrides all other cases
    lob <- plotRange[[varname]][1]
    upb <- plotRange[[varname]][2]
    testlog <- T ## assumes user gave plotRange in canonical scale
  } else if (varname %in% fittedNames) {
    FONKgLow <- blackbox.getOption("FONKgLow")
    FONKgUp <- blackbox.getOption("FONKgUp")
    lob <- FONKgLow[varname]
    upb <- FONKgUp[varname]
    testlog <- F ## FONKgLow is already in krigScale, no need to transform it again
  } else if ( ! is.null(varnameS) && varname %in% varnameS ) { ## to construct a range for a non-Kriging variable but from the Kriged points only
    gridfnchull <- providefullhull(varnameS)[[1]]$vertices ##
    lob <- min(gridfnchull[, varname], na.rm=T)
    upb <- max(gridfnchull[, varname], na.rm=T)
    testlog <- F ## hull was constructed using islogscale: variable is already logscale if islogscale() is true
  } else if (varname %in% blackbox.getOption("ParameterNames")) {
    pointls <- blackbox.getOption("pointls")
    lob <- min(pointls[, varname])
    upb <- max(pointls[, varname])
    testlog <- T ## pointls is in canonical scale
  } else if (varname=="latt2Ns2") {
    latt2Ns2 <- blackbox.getOption("rosglobal")$latt2Ns2
    latt2Ns2pt <- blackbox.getOption("latt2Ns2pt")
    lob <- max(min(latt2Ns2pt), latt2Ns2/2.5)
    upb <- min(max(latt2Ns2pt), latt2Ns2 * 1000)
    testlog <- T
  } else if(varname=="Nratio") {
    Nratiopt <- blackbox.getOption("Nratiopt")
    lob <- min(Nratiopt)
    upb <- max(Nratiopt)
    testlog <- T ## Nratiopt is in canonical scale...
  }  else if(varname=="Nancratio") {
    Nratiopt <- blackbox.getOption("Nratiopt")
    lob <- min(Nratiopt)
    upb <- max(Nratiopt)
    testlog <- T ## Nancratiopt is in canonical scale...
  }  else if(varname=="NactNfounderratio") {
    NactNfounderratiopt <- blackbox.getOption("NactNfounderratiopt")
    lob <- min(NactNfounderratiopt)
    upb <- max(NactNfounderratiopt)
    testlog <- T ## NactNfounderratiopt is in canonical scale...
  }  else if(varname=="NfounderNancratio") {
    NfounderNancratiopt <- blackbox.getOption("NfounderNancratiopt")
    lob <- min(NfounderNancratiopt)
    upb <- max(NfounderNancratiopt)
    testlog <- T ## NfounderNancratiopt is in canonical scale...
  }  else {
    stop.redef(paste("(!) From gridfn: Unhandled variable name '", varname, "'"))
  }
  ##
  if(upb==lob) return(lob) ## single value => to be used as indicator of no variation
  if (testlog && islogscale(varname)) { ## testlog is false if lob, upb values were taken from kriged values
    lob <- log(lob)
    upb <- log(upb)
  }
  if(margefrac>=0.5) message.redef("(!) From gridfn(): 'margefrac' argument too large")
  marge <- margefrac*(upb-lob)
  Grid <- seq(lob+marge, upb-marge, length.out=gridsteps) ## uniform on log scale if relevant
  return(Grid)
}
