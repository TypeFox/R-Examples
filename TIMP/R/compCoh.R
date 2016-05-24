"compCoh" <-
function (irfpar, x, cohspec, coh, dataset, cohirf, mirf = FALSE, 
measured_irf = vector(), convalg = 1, shiftmea = vector(), lamb = 1, 
ani = list(), anipar = vector(), cohcol = vector()) 
{
  type <- cohspec$type 
  if(tolower(type) == "irf") {
    if(mirf) {
      if (length(shiftmea) != 0) {
        if (length(shiftmea) == 1) 
          lamb <- 1
        xspace <- x[2] - x[1]
        measured_irf <- .C("ShiftCurve", 
			   source = as.double(measured_irf), 
			   as.double(measured_irf), 
			   as.double(shiftmea[lamb]/xspace), 
			   as.integer(length(x)), PACKAGE="TIMP")$source
      }
      cohcols <- measured_irf
    }
    else {
      # todo: implement code that takes into account non-gaussian IRF
      cohcols <- dnorm(x, irfpar[1], irfpar[2])
      }
  }
  if(tolower(type) == "freeirfdisp") 
    cohcols <- dnorm(x, cohirf[1], cohirf[2])
  if(tolower(type) == "irfmulti") { 
    cohcols <- matrix(0, nrow = length(x), 
                      ncol = cohspec$numdatasets)
    cohcols[, dataset] <- dnorm(x, irfpar[1], irfpar[2])
  }
  if(tolower(type) == "seq")
    #TODO: fix IRF parameters
    cohcols <- calcCirf(coh, x, irfpar) %*% calcB(coh) 
  if(tolower(type) == "mix") {
    #TODO: fix IRF parameters
    cohcols <- cbind(dnorm(x, irfpar[1], 
			   irfpar[2]), calcCirf(coh, x, irfpar) %*% calcB(coh)) 
  } 
  if(tolower(type) == "xpm") {
    # A*(t*Exp(-2(t/tau)^2))-(t-Tgvd)(Exp(-2((t-Tgvd)^2)/tau^2)
    # A = coh[1]; Tgvd=coh[2]; tau=coh[3] 
    t0 = irfpar[1]
    tt=x-t0
    A=coh[1]
    Tgvd=coh[2]
    tau=coh[3]
    cohcols <- (A*(tt*exp(-2*(tt/tau)^2))-(tt-Tgvd)*(exp(-2*((tt-Tgvd)/tau)^2)))
  }
                                        # constant anisotropy
  if(ani$calcani) {
    if(ani$rammanest) {
      if(mirf) {
        if (length(shiftmea) != 0) {
          if (length(shiftmea) == 1) 
            lamb <- 1
          xspace <- x[2] - x[1]
          measured_irf <- .C("ShiftCurve", 
                             source = as.double(measured_irf), 
                             as.double(measured_irf), 
                             as.double(shiftmea[lamb]/xspace), 
                             as.integer(length(x)), PACKAGE="TIMP")$source
			 }
        ramman <- measured_irf
      }
      else
        ramman <- dnorm(x, irfpar[1], irfpar[2])
      cohcols <- cbind(cohcols, ramman)
	       }
    if(ani$angle[dataset]  != "MA") {
      if(ani$angle[dataset] == "PAR") 
        gamma <- 2
      if(ani$angle[dataset] == "PERP")
        gamma <- -1 
			for(i in 1:length(cohcol))
                          cohcols <- cohcols *
                            (1+(gamma * anipar[ani$parperm[[i]][1]])) 
    }
  }
  cohcols
}

