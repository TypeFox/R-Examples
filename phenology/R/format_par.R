#' .format_par calculates simplified set of parameters.
#' @title The function ".format_par" is for internal use use only.
#' @author Marc Girondot
#' @return Return a modified set of parameters
#' @param xpar Set of parameters
#' @param serie Name of the series to be analyzed
#' @description Calculate a simplified set of parameters.



.format_par <- function(xpar, serie) {
  
#  xpar <- c(Min=12, Peak_Alpha=15, Peak_Beta=-16, Theta=16, Begin=15);serie="Alpha"

  xpar <- na.omit(xpar)
  
#  xpar_courant <<- xpar
#  save.image("courant.RData")

  nxparec <- strsplit(names(xpar), "_")
  ec <- sapply(nxparec, function(x) ifelse(length(x)>1, x[[2]]==serie, TRUE))
  
  xparec <- xpar[ec]
  names(xparec) <- sapply(nxparec[ec], function(x) x[[1]])

if (is.na(xparec["MinB"]) && is.na(xparec["PMinB"]) && is.na(xparec["Min"]) && is.na(xparec["PMin"])) {xparec["MinB"]<-0}
if (is.na(xparec["MinE"]) && is.na(xparec["PMinE"]) && is.na(xparec["Min"]) && is.na(xparec["PMin"])) {xparec["MinE"]<-0}
  
  ec <- !is.na(match(names(xparec), 
                     c("MinB", "MinE", "Min", "PMin", "PminE", "PMinB", "Peak", 
                       "Begin", "End", "Max", "Theta", "Length", "LengthB", "LengthE", 
                       "Tau", "Tau1", "Tau2")))
  xparec[ec] <- abs(xparec[ec])

xparec["Flat"] <- ifelse(is.na(xparec["Flat"]), 0, abs(xparec["Flat"]))

if (!is.na(xparec["Length"])) {
	xparec["Begin"]=xparec["Peak"]-xparec["Length"]
	xparec["End"]=xparec["Peak"]+xparec["Length"]	
	}
if (!is.na(xparec["LengthB"])) {
	xparec["Begin"]=xparec["Peak"]-xparec["LengthB"]
	}
if (!is.na(xparec["LengthE"])) {
	xparec["End"]=xparec["Peak"]+xparec["LengthE"]	
	}



if (!is.na(xparec["PMinE"])) {xparec["MinE"]<-xparec["Max"]*xparec["PMinE"]/100}
if (!is.na(xparec["PMinB"])) {xparec["MinB"]<-xparec["Max"]*xparec["PMinB"]/100}

if (!is.na(xparec["Min"])) {
	xparec["MinB"]<-xparec["Min"]
	xparec["MinE"]<-xparec["Min"]
	}
	
if (!is.na(xparec["PMin"])) {
	xparec["MinB"]<-xparec["Max"]*xparec["PMin"]/100
	xparec["MinE"]<-xparec["Max"]*xparec["PMin"]/100	
	}

xparec["sin"]<-(!is.na(xpar["Phi"]) && !is.na(xpar["Delta"]))
if (xparec["sin"]) {
	if (is.na(xparec["Alpha"])) {xparec["Alpha"]=0}
	if (is.na(xparec["Beta"])) {xparec["Beta"]=0}
	if (is.na(xparec["Tau"])) {xparec["Tau"]=1}
}

xparec["sin1"]<-(!is.na(xpar["Phi1"]) && !is.na(xpar["Delta1"]))
if (xparec["sin1"]) {
	if (is.na(xparec["Alpha1"])) {xparec["Alpha1"]=0}
	if (is.na(xparec["Beta1"])) {xparec["Beta1"]=0}
	if (is.na(xparec["Tau1"])) {xparec["Tau1"]=1}
}

xparec["sin2"]<-(!is.na(xpar["Phi2"]) && !is.na(xpar["Delta2"]))
if (xparec["sin2"]) {
	if (is.na(xparec["Alpha2"])) {xparec["Alpha2"]=0}
	if (is.na(xparec["Beta2"])) {xparec["Beta2"]=0}
	if (is.na(xparec["Tau2"])) {xparec["Tau2"]=1}
}


xparec["PmoinsF"]<-xparec["Peak"]-(xparec["Flat"]/2)
xparec["PplusF"]<-xparec["Peak"]+(xparec["Flat"]/2)

xparec["PmoinsFB"]<-xparec["PmoinsF"]-xparec["Begin"]
xparec["EPplusF"]<-xparec["End"]-xparec["PplusF"]

xparec["MaxMinB"]<-xparec["Max"]-xparec["MinB"]
xparec["MaxMinE"]<-xparec["Max"]-xparec["MinE"]

return(xparec)

}
