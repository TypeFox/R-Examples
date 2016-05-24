calc1DCIs <- function(oneDimCIvars,
                      FONKgNames,
                      fittedNames,
                      CIlevel=blackbox.getOption("CIlevel"),
                      nextBounds=blackbox.getOption("nextBounds"),
                      NextBoundsLevel=blackbox.getOption("NextBoundsLevel"),
                      boundsOutfile="",
                      dataString="",
                      cleanResu=""
) {
  paramnbr <- blackbox.getOption("paramnbr") ## length(ParameterNames)
  if (is.null(paramnbr)) paramnbr <- length(FONKgNames)
  if (is.null(NextBoundsLevel)) NextBoundsLevel <- 0.001
  NextBoundsdlr <- -qchisq(1-NextBoundsLevel, 1)/2
  if (is.null(nextBounds)) nextBounds <- "pointsFromR"
  if(tolower(nextBounds) %in% "frompoints") message.redef("Computing bounds from points...")
  if(tolower(nextBounds) %in% "from1dci") message.redef("Computing bounds from one dim CI...")
  CIlo <- NA;CIup <- NA
  boundstable <- as.data.frame(array(1, c(paramnbr, 4))) ## 1 by default -> values not used
  rosglobal <- blackbox.getOption("rosglobal")
  Nbfactor <- blackbox.getOption("Nbfactor") ## may be NULL
  ## If next migraine iteration depends on from1dci we need to complement oneDimCIvars with any missing variable:
  if (tolower(nextBounds) %in% c("from1dci", "frompoints")) oneDimCIvars <- unique(c(oneDimCIvars, fittedNames))
  if (length(oneDimCIvars)>0) {
    write("\n*** Confidence intervals *** \n", file=cleanResu)
    message.redef("\n*** Confidence intervals *** ")
  }
  CIpointsList <- list()
  for (CIvar in oneDimCIvars) {
    CIvarOri <- CIvar
    if (CIvar=="g" && "condS2" %in% fittedNames) {
      message.redef("");
      message.redef("Requested CI for 'g' while kriging was done for condS2. A CI for condS2 is performed and bounds are converted to 'g' values")
      condS2forgbool <- TRUE
      D2bool <- ("2D" %in% blackbox.getOption("DemographicModel"))
      CIvar <- "condS2"
    } else condS2forgbool <- FALSE
    if ( ("Migraine" %in% blackbox.getOption("usedBy")) && CIvar=="Nb") intlCIvar <- "latt2Ns2" else intlCIvar <- CIvar ## intlCIvar contains the standard internal variable, vs user-style named CIvar
    ###### For each CIvar, optionally outputs CI matching bounds in boundstable if the latter is computed
    ###### Then computes CI according to dedicated CI options
    if (intlCIvar %in% FONKgNames) rownames(boundstable)[which(intlCIvar==FONKgNames)] <- intlCIvar
    ## here we seek values that are either directly written in boundstable ('frompoints' case)
    ## or are starting points and bounds for profilind1Dwrapper
    if( ! (intlCIvar %in% fittedNames)) {
      if(intlCIvar=="latt2Ns2") {
        MLval <- NA ## should be recomputed using composite hull within bounds1D
      } else if(CIvar=="Nratio") {
        MLval <- NA ## should be recomputed using composite hull within bounds1D
      } else if(CIvar=="Nancratio") {
        MLval <- NA ## should be recomputed using composite hull within bounds1D
      } else if(CIvar=="NactNfounderratio") {
        MLval <- NA ## should be recomputed using composite hull within bounds1D
      } else if(CIvar=="NfounderNancratio") {
        MLval <- NA ## should be recomputed using composite hull within bounds1D
      } else if(CIvar=="twoNm") {
        MLval <- NA ## should be recomputed using composite hull within bounds1D
      } else {
        immedst <- paste("(!) 1D CI computation on given combination of CI/kriging variables not implemented : ", CIvar)
        message.redef(immedst)
        write(immedst, file=cleanResu)
        ## ... in particular no providefulhull algo in this case...
        civarst <- paste(dataString, "(", CIvar, "_CI_NOT_COMPUTED)", sep="")
        writeoutput(civarst=civarst, levelSlot=NA, CIloSlot=NA, CIupSlot=NA)
        next ## SKIPS remainder of loop body
      }
    } else { ## CI variable in kriging space
      ## first composite kriging var
      if (intlCIvar=="latt2Ns2") {
        MLval <- rosglobal$latt2Ns2
      } else if (intlCIvar=="Nratio") {
        MLval <- rosglobal$Nratio ## FR: written by analogy, maybe not meaningful
      } else if (intlCIvar=="Nancratio") {
        MLval <- rosglobal$Nratio ## RL: written by analogy, maybe not meaningful
      } else if (intlCIvar=="NactNfounderratio") {
        MLval <- rosglobal$NactNfounderratio ## RL: written by analogy, maybe not meaningful
      } else if (intlCIvar=="NfounderNancratio") {
        MLval <- rosglobal$NfounderNancratio ## RL: written by analogy, maybe not meaningful
      } else { ## canonical kriging var
        MLval <- (rosglobal$canonVP)[intlCIvar]
      }
      if (islogscale(intlCIvar)) MLval <- log(MLval)
    }
    locchull.pts <- providefullhull(intlCIvar)[[1]]$vertices
    if (! intlCIvar %in% colnames(locchull.pts)) {
      stop.redef("Unfeasible CI requested; either dubious combination of options (e.g., CI for 'Nb' requested with 'g' fixed) or a bug.")
    }
    lowval <- min(locchull.pts[, intlCIvar])
    hival <- max(locchull.pts[, intlCIvar])
    ## bloc pour boundstable
    if (tolower(nextBounds) %in% c("from1dci", "frompoints")) {
      if(tolower(nextBounds)=="frompoints") {
        CIlo <- lowval
        CIup <- hival
      } else if(tolower(nextBounds)=="from1dci") { ##
        CIloup <- bounds1D(NextBoundsdlr, "Next lower bound from 1DCI", "Next upper bound from 1DCI", lowval, hival, MLval, intlCIvar)
        CIlo <- CIloup$CIlo;CIup <- CIloup$CIup
      }
      if (islogscale(intlCIvar)) {
        CIlo <- exp(CIlo)
        CIup <- exp(CIup)
      }
      boundstable[intlCIvar, ] <- c(is.na(CIlo), is.na(CIup), ifelse(is.na(CIlo), 666, CIlo), ifelse(is.na(CIup), 666, CIup))
      ## write NextBoundsLevel CI to output[...].txt file too! (perhaps not a good idea with frompoints)
      if (CIvar=="Nb") {
        if (! is.null(Nbfactor)) {
          CIlo <- CIlo*Nbfactor
          CIup <- CIup*Nbfactor
        }
      } ## would be the standard user's practice
      if (condS2forgbool) {
        CIlo <- groot(CIlo, D2bool=D2bool)
        CIup <- groot(CIup, D2bool=D2bool)
      }
      civarst <- paste(dataString, "(", CIvarOri, "_CI)", sep="")
      writeoutput(civarst=civarst, levelSlot=NextBoundsLevel, CIloSlot=CIlo, CIupSlot=CIup)
    }
    ## bloc not for boundstable
    ##    for 1DCI with CIdlr (not nextBoundsdlr, so by def does not go in nextbounds)
    CIdlr <- -qchisq(1-CIlevel, 1)/2
    if ("NLOPT_LN_COBYLA_for_CI" %in% blackbox.getOption("optimizers")
        || "NLOPT_LD_MMA_for_CI" %in% blackbox.getOption("optimizers")) {
      CIloup <- calcBounds1D(CIdlr, "Lower CI bound", "Upper CI bound", lowval, hival, maxparval=MLval, intlCIvar)
    } else {
      CIloup <- bounds1D(CIdlr, "Lower CI bound", "Upper CI bound", lowval, hival, MLval, intlCIvar)
    }
    CIlo <- CIloup$CIlo
    CIup <- CIloup$CIup
    CIpointsList[[intlCIvar]] <- CIloup$CIpoints ## je retire les noms comme ça les elements n'heriterons que le nom de la variable
    if (islogscale(intlCIvar)) {CIlo <- exp(CIlo);CIup <- exp(CIup)}
    if (CIvar=="Nb") {
      if (! is.null(Nbfactor)) {
        CIlo <- CIlo*Nbfactor
        CIup <- CIup*Nbfactor
      }
    } ## would be the standard user's practice
    if (condS2forgbool) {CIlo <- groot(CIlo, D2bool=D2bool);CIup <- groot(CIup, D2bool=D2bool)}
    do.call(blackbox.options, list(CIlo= CIlo, CIup=CIup))
    username <- userunit(CIvarOri, format="ASCII")
    civarst <- paste(dataString, "(", username, "_CI)", sep="")
    writeoutput(civarst=civarst, levelSlot=CIlevel, CIloSlot=CIlo, CIupSlot=CIup)
    immedst <- paste(100*(1-CIlevel), "%-coverage confidence interval for ",
                     username, " : [ ", prettynum(CIlo), " -- ", prettynum(CIup), " ]", sep="")
    write(immedst, file=cleanResu)
  } ## end for (CIvar in oneDimCIvars)
  if (tolower(nextBounds) %in% c("from1dci", "frompoints")) {
    write(t(boundstable), file=boundsOutfile, ncolumns=4)
  }
  return(CIpointsList)
}
