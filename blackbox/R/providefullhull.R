setHullPrecMode <- function(notinKgspace,outputnames) {
  namebit <- ""
  optRedmode <- blackbox.getOption("redundant.mode")
  if ("latt2Ns2" %in% notinKgspace) {
    outputnames[which(outputnames=="twoNm")] <- "latt2Ns2"
    namebit <- paste(namebit, "Ns2", sep="_")
    if (optRedmode=="defaultPrecision") redmode <- "rational"
  }
  if ("twoNm" %in% notinKgspace) {
    outputnames[which(outputnames=="latt2Ns2")] <- "twoNm"
    namebit <- paste(namebit, "Nm", sep="_")
    if (optRedmode=="defaultPrecision") redmode <- "double"
  }
  if ("condS2" %in% notinKgspace) {
    outputnames[which(outputnames=="g")] <- "condS2"
    namebit <- paste(namebit, "condS2", sep="_")
    if (optRedmode=="defaultPrecision") redmode <- "double"
  }
  if ("g" %in% notinKgspace) {
    outputnames[which(outputnames=="condS2")] <- "g"
    namebit <- "NEVER_CONSIDERED_YET"
    if (optRedmode=="defaultPrecision") redmode <- "double"
  }
  if ("Nratio" %in% notinKgspace) {
    outputnames[which(outputnames=="twoNmu")] <- "Nratio"
    namebit <- paste(namebit, "Nratio", sep="_")
    if (optRedmode=="defaultPrecision") redmode <- "rational"
  }
  if ("NactNfounderratio" %in% notinKgspace) {
    outputnames[which(outputnames=="twoNmu")] <- "NactNfounderratio"
    namebit <- paste(namebit, "NactNfounderratio", sep="_")
    if (optRedmode=="defaultPrecision") redmode <- "rational"
  }
  if ("NfounderNancratio" %in% notinKgspace) {
    outputnames[which(outputnames=="twoNancmu")] <- "NfounderNancratio"
    namebit <- paste(namebit, "NfounderNancratio", sep="_")
    if (optRedmode=="defaultPrecision") redmode <- "rational"
  }
  if ("Nancratio" %in% notinKgspace) {
    outputnames[which(outputnames=="twoNmu")] <- "Nratio"
    namebit <- paste(namebit, "Nratio", sep="_")
    if (optRedmode=="defaultPrecision") redmode <- "rational"
  }
  if(nchar(namebit)==0L) {stop.redef("composite variable not handled in hull computation")}
  return(list(redmode=redmode,outputnames=outputnames,namebit=namebit))
}



## pointsinKgSpace provides for computation of hull from given point.
## If it is used, no input/output to options$hulls should be used
# The pointsinKgSpace variable must not be changed within the code since the original nullness is tested
#
providefullhull <- function(varnames) { ##varnames should include variables not in FONKgNames but may also contain variables in FONKgNames
  fittedNames <- blackbox.getOption("fittedNames")
  redmode <- switch(blackbox.getOption("redundant.mode"),
           "noElim"="no.elim",
           "alwaysRational"="rational",
           "alwaysDouble"="double",
           "defaultPrecision"= if (length(blackbox.getOption("varnames") %w/o%
                                          blackbox.getOption("ParameterNames"))>0) {"rational"} else {"double"}
    )
  ## default method depends on notinKgspace below
  notinKgspace <- varnames %w/o% blackbox.getOption("FONKgNames")
  tmp1 <- blackbox.getOption("FONKgpointls")[, fittedNames, drop=FALSE]
  if (length(notinKgspace)==0) {
    if (is.null(locchull <- blackbox.getOption("hulls")$Kgtotal)) {
      if (  blackbox.getOption("redundant.mode")=="defaultPrecision") redmode <- "double"
      locchull <- resetCHull(tmp1, formats=c("vertices", "vertices001", "constraints"), redundant.mode=redmode)
      .blackbox.data$options$hulls$Kgtotal <- locchull ## ici l'acces direct aux membres de la liste est utile
    }
    return(list(Kgtotal=locchull))
  }
  ##ELSE
  ## verif validity of varnames and builds outputnames
  outputnames <- fittedNames
  if (length(notinKgspace)>0) {
    locblob <- setHullPrecMode(notinKgspace=notinKgspace,
                               outputnames=outputnames)
    outputnames <- locblob$outputnames
    redmode <- locblob$redmode
    namebit <- locblob$namebit
  }
  locchull <- blackbox.getOption("hulls")[[namebit]] ## may be NULL
  if (is.null(locchull)) {
    tmp1 <- t(apply(tmp1, 1, fromFONKtoanyspace, outputnames=outputnames))
    colnames(tmp1) <- outputnames  ## the outputname are lost by apply which keeps the original names !!
    locchull <- resetCHull(tmp1, formats=c("vertices", "constraints"), redundant.mode=redmode)
    locchull <- matchVertCons(locchull) ## ADDS correspondance between vertices and constraints
    .blackbox.data$options$hulls[[namebit]] <- locchull ## .blackbox required here
  }
  resu <- list(locchull)
  names(resu) <- namebit
  return(resu)
}


provideVertices <- function(varnames, pointsinKgSpace) {
  fittedNames <- blackbox.getOption("fittedNames")
  redmode <- switch(blackbox.getOption("redundant.mode"),
           "noElim"="no.elim",
           "alwaysRational"="rational",
           "alwaysDouble"="double",
           "defaultPrecision"= if (length(blackbox.getOption("varnames") %w/o%
                                          blackbox.getOption("ParameterNames"))>0) {"rational"} else {"double"}
    )
  notinKgspace <- varnames %w/o% blackbox.getOption("FONKgNames")
  ## verif validity of varnames and builds outputnames
  outputnames <- fittedNames
  if (length(notinKgspace)>0) {
    locblob <- setHullPrecMode(notinKgspace=notinKgspace,
                               outputnames=outputnames)
    outputnames <- locblob$outputnames
    redmode <- locblob$redmode
  }
  inoutspace <- apply(pointsinKgSpace, 1, fromFONKtoanyspace, outputnames=outputnames)
  if (length(outputnames)>1L) {
    inoutspace <- t(inoutspace)
  } else inoutspace <- matrix(inoutspace, ncol=1)
  colnames(inoutspace) <- outputnames  ## the outputname are lost by apply which keeps the original names !!
  colmins <- apply(inoutspace, 2, min)
  colmaxs <- apply(inoutspace, 2, max)
  ##Note that FONKgLow/Up will be recomputed one the points have been selected for Kriging
  fixBools <- ((colmaxs-colmins)<1e-06) ##FR->FR test pas compar a range de FONKg ?
  if (any(fixBools)) inoutspace <- inoutspace[, ( ! fixBools), drop=FALSE] ## do not provide constant cols to resetCHull->convhulln...
  locvertices <- resetCHull(inoutspace, formats=c("vertices"), redundant.mode=redmode)$vertices
  resu <- matrix(NA,ncol=length(outputnames),nrow=nrow(locvertices))
  colnames(resu) <- outputnames
  resu[,colnames(inoutspace)] <- locvertices
  for (st in outputnames[fixBools]) resu[,st] <- colmins[st]
  return(resu)
}
