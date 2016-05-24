userunit <- function(charname, logscalest="", value, format="expression", Nbstyle=T) { ##
  tmp <- formatName(charname, format=format, Nbstyle=Nbstyle)
  fullstring <- logscalest
  if ( ("1D" %in% blackbox.getOption("DemographicModel")) && charname=="latt2Ns2") {
    if (Nbstyle) {
      fullstring <- paste(" (", blackbox.getOption("GeoUnit"), ")", fullstring, sep="") ## default for Nb in 1D
    } else fullstring <- paste(" (lattice units)", fullstring, sep="") ##
  }
  if (format=="expression") fullstring <- chartr(" ", "~~", fullstring)
  if(format!="expression") {return(paste(tmp, fullstring, sep=""))} else {return(parse(text=paste(tmp, eval(bquote(.(fullstring))), sep="")))}
}
