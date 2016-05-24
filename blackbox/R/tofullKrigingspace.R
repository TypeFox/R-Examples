tofullKrigingspace <- function(fittedlist, fixedlist=NULL) { ## output in kriging space (reorder columns in standard order if required)
  fittedparamnbr <- blackbox.getOption("fittedparamnbr")
  fittedNames <- blackbox.getOption("fittedNames")
  ParameterNames <- blackbox.getOption("ParameterNames")
  getvalue <- function(st) {
    if (st %in% names(fittedlist)) {
      return(as.numeric(fittedlist[st]))
    } else if (st %in% names(fixedlist)) {
      return(as.numeric(fixedlist[st]))
    } else if (st %in% blackbox.getOption("constantNames")) { ## we may need values of fixed parameters
      return(as.numeric(blackbox.getOption("FONKgLow")[st]))
    } else return(NA)
  }
  ## fittedlist and fixedlist must be in FONKgScale or else in extraScale
  ## (01/2010: 'full' Krig sp of *fittedparamnbr* (the argument of purefn)
  ## Generates a 'full' vector from fittedlist and fixedlist. fixedlist can be NULL although few subcases are implemented
  ## check fittedlist
  if(is.null(names(fittedlist))) {
    stop.redef("(!) names(fittedlist) is NULL in tofullKrigingspace. ")
  }
  ## FR 18/11/10 replaced fittedNames by ParameterNames in next line, for case where latt2Ns2 in fittedNames and twoNm in fixedlist (for LRT on twoNm...)
  ## Not sure whether this will be OK with fittedparamnbr< param nbr
  checkfitted <- names(fittedlist) %w/o% c(ParameterNames, "latt2Ns2", "Nratio", "NactNfounderratio", "NfounderNancratio", "condS2")
  if (length(checkfitted)>0) {
    message.redef("(!)From tofullKrigingspace(): names(fittedlist) contains unhandled variable")
    message.redef("   or combination of variables") ##latt2Ns2 plus another variable
    message.redef(checkfitted) ## hind: dont forget to use <name>=<vector>[[<name>]], ie [[ ]] not [ ]
    stop.redef()
  }
  ## check fixedlist
  if(!is.null(fixedlist)) {
    checkfixed <- names(fixedlist) %w/o% c(ParameterNames, "latt2Ns2", "Nratio", "NactNfounderratio", "NfounderNancratio", "condS2")
    if (length(checkfixed)>0) {
      message.redef("(!)From tofullKrigingspace(): names(fixedlist) contains unhandled variable")
      message.redef("   or combination of variables") ##latt2Ns2 plus another variable
      message.redef(paste(names(fixedlist)))
      stop.redef()
    }
  }
  KrigVec <- rep(NA, fittedparamnbr)
  names(KrigVec) <- fittedNames
  ## if everything is in canonical space then the next two lines fill the KrigVec vector
  finf <- intersect(names(fixedlist), fittedNames)
  KrigVec[finf] <- as.numeric(fixedlist[finf])
  finf <- intersect(names(fittedlist), fittedNames)
  KrigVec[finf] <- as.numeric(fittedlist[finf])
  if (!any(is.na(KrigVec))) { ## then we are done
    return(KrigVec)
  } ##ELSE
  ## we have to handle composite variables. We first unlog everything that is logscale
  ## we must be careful that even if say extrascale=Nb=logscale, a non-log latt2Ns2 can be given in the arguments to tofullKrigingspace...
  for(st in names(fittedlist)) if (islogscale(st)) {fittedlist[[st]] <- exp(fittedlist[[st]])}
  for(st in names(fixedlist)) if (islogscale(st)) {fixedlist[[st]] <- exp(fixedlist[[st]])}
  for(st in names(KrigVec)) if (islogscale(st)) {KrigVec[[st]] <- exp(KrigVec[[st]])} ## exp(NA)->NA, not a problem
  D2bool <- ("2D" %in% blackbox.getOption("DemographicModel"))
  ## then we operate in canonical scale
  if("twoNm" %in% fittedNames && is.na(KrigVec["twoNm"])) { ##
    latt2Ns2 <- getvalue("latt2Ns2")
    if(is.na(latt2Ns2)) stop.redef("(!) From tofullKrigingspace(): neither twoNm nor latt2Ns2 given")
    g <- getvalue("g")
    if (is.na(g)) {
      S2 <- getvalue("condS2")
      if (is.na(S2)) stop.redef("(!) From tofullKrigingspace(): neither twoNm nor g nor condS2 given")
    } else S2 <- condaxialS2fromg(g, D2bool=D2bool)
    KrigVec["twoNm"] <- latt2Ns2/S2
  }
  if("latt2Ns2" %in% fittedNames && is.na(KrigVec["latt2Ns2"])) { ## then we must have twoNm somewhere
    twoNm <- getvalue("twoNm")
    if(is.na(twoNm)) stop.redef("(!) From tofullKrigingspace(): neither twoNm nor latt2Ns2 given")
    S2 <- getvalue("condS2")
    if (is.na(S2)) {
      S2 <- condaxialS2fromg(KrigVec["g"], D2bool=D2bool) ##cond axial S2 from g
    }
    KrigVec["latt2Ns2"] <- twoNm*S2 ## 2Nm as latt2Ns2/condS2
  }
  if("condS2" %in% fittedNames && is.na(KrigVec["condS2"])) {
    g <- getvalue("g")
    if(is.na(g)) { ## need to reconstruct it from other information
      twoNm <- getvalue("twoNm")
      if(is.na(twoNm)) stop.redef("(!) From tofullKrigingspace(): neither twoNm nor condS2 given")
      latt2Ns2 <- getvalue("latt2Ns2")
      if(is.na(latt2Ns2))  stop.redef("(!) From tofullKrigingspace(): neither twoNm nor latt2Ns2 given")
      KrigVec["condS2"] <- latt2Ns2/twoNm
    } else KrigVec["condS2"] <- condaxialS2fromg(g, D2bool=D2bool) ##cond axial S2 from g
  }
  if("g" %in% fittedNames && is.na(KrigVec["g"])) { ## then we must have twoNm and latt2Ns2 somewhere
    twoNm <- getvalue("twoNm")
    S2 <- getvalue("condS2")
    if (is.na(S2)) {
      latt2Ns2 <- getvalue("latt2Ns2")
      if (is.na(twoNm) || is.na(latt2Ns2) ) stop.redef("(!) From tofullKrigingspace(): g and [either twoNm or latt2Ns2] are not given")
      S2 <- latt2Ns2/twoNm ##cond axial S2 from g
    }
    KrigVec["g"] <- groot(S2, D2bool=D2bool )
  }
  if("twoNmu" %in% fittedNames && is.na(KrigVec["twoNmu"])) { ## then we must have Nratio somewhere
    Nratio <- getvalue("Nratio")
    if(is.na(Nratio)) {
      NactNfounderratio <- getvalue("NactNfounderratio")
      if(is.na(NactNfounderratio)) stop.redef("(!) From tofullKrigingspace(): neither twoNmu nor Nratio nor NactNfounderratio given")
      if(is.na(KrigVec["twoNfoundermu"])) { ## RL 052013 means that twoNfoundermu was not fitted
        KrigVec["twoNmu"] <- NactNfounderratio*blackbox.getOption("FONKgLow")["twoNfoundermu"]
      } else {KrigVec["twoNmu"] <- NactNfounderratio*KrigVec["twoNfoundermu"]} ## twoNfoundermu is already unlog'ed above and twoNmu will be relog'ed below
    } else {
      if(is.na(KrigVec["twoNancmu"])) { ## RL 052013 means that twoNancmu was not fitted
        KrigVec["twoNmu"] <- Nratio*blackbox.getOption("FONKgLow")["twoNancmu"]
      } else {KrigVec["twoNmu"] <- Nratio*KrigVec["twoNancmu"]} ## twoNancmu is already unlog'ed above and twoNmu will be relog'ed below
    }
  }
  if("twoNancmu" %in% fittedNames && is.na(KrigVec["twoNancmu"])) { ## then we must have Nratio somewhere
    Nratio <- getvalue("Nratio")
    if(is.na(Nratio)) {
      NfounderNancratio <- getvalue("NfounderNancratio")
      if(is.na(NfounderNancratio)) stop.redef("(!) From tofullKrigingspace(): neither twoNancmu nor Nratio nor NfounderNancratio given")
      if(is.na(KrigVec["twoNfoundermu"])) { ## RL 022016 means that twoNfoundermu was not fitted
        KrigVec["twoNancmu"] <- blackbox.getOption("FONKgLow")["twoNfoundermu"]/NfounderNancratio
      } else {KrigVec["twoNancmu"] <- KrigVec["twoNfoundermu"]/NfounderNancratio} ## twoNfoundermu is already unlog'ed above and twoNmu will be relog'ed below
    } else {
      if(is.na(KrigVec["twoNmu"])) { ## RL 022016 means that twoNmu was not fitted
        KrigVec["twoNancmu"] <- blackbox.getOption("FONKgLow")["twoNmu"]/Nratio
      } else {KrigVec["twoNancmu"] <- KrigVec["twoNmu"]/Nratio} ## twoNmu is already unlog'ed above and twoNancmu will be relog'ed below
    }
  }
  if(!all(is.numeric(KrigVec))) {
    message.redef("(!) From tofullKrigingspacefn(): !all(is.numeric(KrigVec))")
    llocalst <- paste("fittedlist was ", fittedlist)
    message.redef(llocalst)
    llocalst <- paste("fixedlist was ", fixedlist)
    message.redef(llocalst)
    llocalst <- paste("KrigVec was ", KrigVec)
    message.redef(llocalst)
  }
  ## Finally we relog everything that is logscale
  for(st in names(KrigVec)) if ((st %in% fittedNames) && islogscale(st)) {KrigVec[[st]] <- log(KrigVec[[st]])}
  return(KrigVec)
} ## end tofullKrigingspace()
