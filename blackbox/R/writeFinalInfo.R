writeFinalInfo <- function(cleanResu="") {
  rosglobal <- blackbox.getOption("rosglobal")
  message.redef("...done.")
  message.redef("\n*** Final likelihood estimates, and predicted logL: ***")
  message.redef(prettyNamedUserValues(c(rosglobal$canonVP, "ln(L)"=rosglobal$value), extradigits=2))
  returncode <- rosglobal$convergence
  tmp <- rosglobal$edgelevel
  if (tmp>0) returncode <- returncode+tmp/(10^ceiling(log(tmp, 10))) ## second summand goes in decimal part of returcode
  write("\n*** Point estimates *** \n", file=cleanResu)
  writeCleanNAMED(prettyNamedUserValues(rosglobal$canonVP), file=cleanResu)
  DemographicModel <- blackbox.getOption("DemographicModel")
  if ("IBD" %in% DemographicModel) {
    D1IBDbool <- "1D" %in% DemographicModel
    Nbfactor <- blackbox.getOption("Nbfactor")
    write(paste("\n      Neighborhood: ", prettynum(rosglobal$latt2Ns2*Nbfactor), " ",
                if (D1IBDbool) {blackbox.getOption("GeoUnit")} else {""}, sep=""), file=cleanResu)
    if (D1IBDbool) write(paste("## Conversion factor from Nb in lattice units to Nb in geographic distance units as deduced from input:", Nbfactor),
                         file=blackbox.getOption("estimOutf"))
  } else if (length(intersect(DemographicModel, c("OnePopVarSize", "IM")))>0) {
    write(paste("\n      N ratio: ", prettynum(rosglobal$Nratio), " ", sep=""), file=cleanResu)
  } else if ("OnePopFounderFlush" %in% DemographicModel) {
    write(paste("\n      Nanc ratio: ", prettynum(rosglobal$Nratio), " ", sep=""), file=cleanResu)
    write(paste("\n      NactNfounder ratio: ", prettynum(rosglobal$NactNfounderratio), " ", sep=""), file=cleanResu)
    write(paste("\n      NfounderNanc ratio: ", prettynum(rosglobal$NfounderNancratio), " ", sep=""), file=cleanResu)
  }
  ## note that the C codes seeks estimates in the VERY LAST line of the output.txt file: do not write comments after the following output:
  ## GOP <- ...RMSpred/...RMSy ## Goodness of prediction criterion.
  improbableError <- blackbox.getOption("RMSpred")*qnorm(0.99, 0, 1)
  GOP <- rosglobal$value + blackbox.getOption("maxobsFONKy") - improbableError
  writeoutput(paste(blackbox.getOption("dataFile"), "(final)", sep=""),
              returncode=returncode,
              prettynum(blackbox.getOption("RMSpred")),
              prettynum(blackbox.getOption("RMSy")),
              prettynum(GOP))
  if ( ! blackbox.getOption("interactiveGraphics")) { ##
    plotFiles <- blackbox.getOption("plotFiles")
    if(!is.null(plotFiles) & length(plotFiles)>0)
      message.redef(paste("See file(s) ", paste(names(plotFiles), sep="", collapse=", "), " for figures", sep=""))
    graphics.off() ## close all graphic files, but note that plotFiles is not cleaned.
    close(blackbox.getOption("estimOutf")) ## ! leaves connection open if interactive run !
  }
  write("\nNormal ending.", file=cleanResu)
  close(cleanResu) ##
  blackbox.options(cleanResu="") ## so that all 'cleanResu' output now goes to the standard output connection; see ?write
  alarm() ## Sounds when execution of the R source file reaches this point. Except in Rstudio...
  invisible(NULL)
}
