canonize <- function(input) { ## from vector within Kriging space AND in Kriging scale, also returns composite var
  FONKinput <- blackbox.getOption("FONKgLow")
  FONKinput[names(input)] <- input
  FONKinput <- unlist(FONKinput) ## previous line sometimes creates a list although no argument is ??
  ## unlogs
  logs <- tolower(.blackbox.data$options$FONKgScale)=="logscale" ## vector of T/F
  FONKinput[logs] <- exp(FONKinput[logs]) ## argument of from2Ns2Tocanon should not be logscale....
  canon <- FONKinput
  names(canon) <- blackbox.getOption("ParameterNames") ## anticipate future names
  DemographicModel <- blackbox.getOption("DemographicModel")
  FONKgNames <- blackbox.getOption("FONKgNames")

  if ("IBD" %in% DemographicModel) {
    if("condS2" %in% FONKgNames) {
      canon["g"] <- groot(FONKinput["condS2"], D2bool= ("2D" %in% DemographicModel) )
    }
    if("latt2Ns2" %in% FONKgNames) {
      latt2Ns2 <- FONKinput["latt2Ns2"] ## saved in return value of canonize
      canon["twoNm"] <- from2Ns2Tocanon(FONKinput)["twoNm"]
    } else {## constructs 2Ds2 [lattice units]
      latt2Ns2 <- (tolatt2Ns2(canon))["latt2Ns2"] ## requires that canon is indeed already canonical
    }
    return(list(canonVP=canon, latt2Ns2=latt2Ns2))
  } else if ( length(intersect(DemographicModel, c("OnePopVarSize", "IM")))>0) {
    if("Nratio" %in% FONKgNames) {
      Nratio <- FONKinput["Nratio"] ## saved in return value of canonize
      canon["twoNmu"] <- FONKinput[["Nratio"]]*FONKinput[["twoNancmu"]]
    } else {## constructs Nratio
      Nratio <- toNratioFromCanonical(canon) ## requires that canon is indeed already canonical
    }
    return(list(canonVP=canon, Nratio=Nratio))
  } else if ("OnePopFounderFlush" %in% DemographicModel) {
    ##Remember that Nratio in FounderFush is called Nancratio for the user, RL 052013
    if("Nratio" %in% FONKgNames) {
      Nratio <- FONKinput["Nratio"] ## saved in return value of canonize
      canon["twoNmu"] <- FONKinput[["Nratio"]]*FONKinput[["twoNancmu"]]
    } else {## constructs Nratio
      Nratio <- toNratioFromCanonical(canon) ## requires that canon is indeed already canonical
    }
    if("NactNfounderratio" %in% FONKgNames) {
      NactNfounderratio <- FONKinput["NactNfounderratio"] ## saved in return value of canonize
      canon["twoNmu"] <- FONKinput[["NactNfounderratio"]]*FONKinput[["twoNfoundermu"]]
    } else {## constructs NactNfounderratio
      NactNfounderratio <- toNactNfounderratioFromCanonical(canon) ## requires that canon is indeed already canonical
    }
    if("NfounderNancratio" %in% FONKgNames) {
      NfounderNancratio <- FONKinput["NfounderNactratio"] ## saved in return value of canonize
      canon["twoNancmu"] <- FONKinput[["twoNfoundermu"]]/FONKinput[["NfounderNancratio"]]
    } else {## constructs NfounderNancratio
      NfounderNancratio <- toNfounderNancratioFromCanonical(canon) ## requires that canon is indeed already canonical
    }
    return(list(canonVP=canon, Nratio=Nratio, NactNfounderratio=NactNfounderratio, NfounderNancratio=NfounderNancratio))
  } else return(list(canonVP=canon))
} ## end def canonize
