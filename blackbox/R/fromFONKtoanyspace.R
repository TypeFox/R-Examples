fromFONKtoanyspace <- function(FONKinput, outputnames, outputScale=blackbox.getOption("FONKgScale")) {
  ## 'FONKinput' is vector in FONKriging space AND SCALE,
  #                          10/2015 FR->FR nom de la fonction bizarre, vu qu'elle convertit from fitted to any space...
  ## 'outputnames' are the return variables
  ## 'outputScale' added 03/2015 to convert to sampling scale when converting to samplingSpace
  ##   its default may see odd, but this was default assumed by the islogcale call without outputScale
  ##   so this retains the older behaviour when there is no explicit outputScale args
  if(is.matrix(FONKinput)) stop.redef("(!) From fromFONKtoanyspace(): this function does not take matrix arguments")
  inputnames <- names(FONKinput)
  fittedNames <- blackbox.getOption("fittedNames")
  if (length(inputnames)!= length(fittedNames) || ! all(inputnames==blackbox.getOption("fittedNames"))) {
    message.redef("(!) From fromFONKtoanyspace(): names of input:")
    message.redef(paste(inputnames,collapse=", "))
    message.redef("    do not match the expected names: ")
    message.redef(paste(fittedNames,collapse=", "))
    stop.redef()
  } ## ELSE
  if (length(outputnames %w/o% inputnames)>0) {
    canon <- canonize(FONKinput) ## unlogs...
    output <- canon$canonVP ## bug: ;names(output) <- outputnames canonVP can be longer than outputnames
    if ("latt2Ns2" %in% outputnames) {names(output)[which(names(output)=="twoNm")] <- "latt2Ns2"; output["latt2Ns2"] <- canon$latt2Ns2}
    if ("condS2" %in% outputnames) {names(output)[which(names(output)=="g")] <- "condS2"; output["condS2"] <- canon$latt2Ns2/canon$canonVP["twoNm"]} ## writes condS2 in place of g
    if ("Nratio" %in% outputnames) {names(output)[which(names(output)=="twoNmu")] <- "Nratio"; output["Nratio"] <- canon$Nratio}
    if ("Nancratio" %in% outputnames) {message.redef("Nancratio Here 7");names(output)[which(names(output)=="twoNmu")] <- "Nratio"; output["Nratio"] <- canon$Nratio}
    if ("NactNfounderratio" %in% outputnames) {names(output)[which(names(output)=="twoNmu")] <- "NactNfounderratio"; output["NactNfounderratio"] <- canon$NactNfounderratio}
    if ("NfounderNancratio" %in% outputnames) {names(output)[which(names(output)=="twoNancmu")] <- "NfounderNancratio"; output["NfounderNancratio"] <- canon$NfounderNancratio}
  } else { ## need to unlog since no call to canonize.
    ## the else close is slightly inelegant (one could rather compare the two logscales), but clear.
    output <- FONKinput
    for(st in inputnames) {if (islogscale(st)) {output[st] <- exp(output[st])}}
  }
  ## at this point the names must match the output names
  ## relogs
  for(st in outputnames) {if (islogscale(st,scale=outputScale)) {output[st] <- log(output[st])}}
  if (any(is.nan(output))) { ### should never occur in normal use, not handled by calling function
    message.redef("(!) From fromFONKtoanyspace(): NaN's in putative return value")
    message.redef(paste("output= ", output))
    stop.redef()
  }
  return(output[outputnames])
} ## end fromFONKtoanyspace()
