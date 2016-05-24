writeoutput <- function(civarst, returncode=NA,
                        levelSlot,
                        CIloSlot,
                        CIupSlot) {
  rosglobal <- blackbox.getOption("rosglobal")
  begstv <- c(civarst, 1, blackbox.getOption("lenptls"), rosglobal$canonVP)
  endstv <- c(rosglobal$value, returncode, levelSlot, CIloSlot, CIupSlot)
  ## plus maxp$par/value instead of rosglobal$par/value for CIvar not in Kg space... ??
  if ("IBD" %in% blackbox.getOption("DemographicModel")) {
    if (is.na(rosglobal$latt2Ns2)) {composite <- NA} else {composite <- rosglobal$latt2Ns2*blackbox.getOption("Nbfactor")}
  } else {
    composite <- ""
  }
  stv <- c(begstv, composite, endstv)
  ncolumns <- blackbox.getOption("paramnbr")+9
  write(stv, file=blackbox.getOption("estimOutf"), ncolumns=ncolumns)
  ##return(canonVPoutput)
} ## end def writeoutput()
