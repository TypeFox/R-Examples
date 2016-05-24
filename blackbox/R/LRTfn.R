LRTfn <- function(LRTfixedvals, cleanResu) { ## LRTfixedvals is already in logscale and has standardNames
  fitobject <- blackbox.getOption("fitobject")
  fittedNames <- blackbox.getOption("fittedNames")
  fittedparamnbr <- blackbox.getOption("fittedparamnbr")
  ### only for outputs to output.txt, and to the screen
  usernames <- sapply(names(LRTfixedvals), userunit, format="ASCII") ## may be nicer than those actually entered by user
  uservalues <- unlist(LRTfixedvals)
  for(st in names(uservalues)) {if (islogscale(st)) {uservalues[st] <- exp(uservalues[st])}}
  uservalues <- prettynum(uservalues, extradigits=10) ## the digits as entered by the users... eg at least two extradigits needed for 0.99999 not shown as 1.
  ###
  locst <- paste(blackbox.getOption("dataFile"), "(LRT_", paste(usernames, uservalues, sep="_", collapse="_"), ")", sep="")
  notinKgspace <- names(LRTfixedvals) %w/o% fittedNames
  if ("g" %in% notinKgspace && "condS2" %in% fittedNames) {
    message.redef("Requested LRT for 'g' while kriging was done for condS2. The LRT for g is deduced from the equivalent LRT for condS2.")
    condS2forgbool <- T
    notinKgspace <- notinKgspace %w/o% "g"
    LRTfixedvals[["condS2"]] <- condaxialS2fromg(LRTfixedvals[["g"]], D2bool=("2D" %in% blackbox.getOption("DemographicModel")))
    if (islogscale("condS2")) LRTfixedvals[["condS2"]] <- log(LRTfixedvals[["condS2"]])
    LRTfixedvals[["g"]] <- NULL
  } else condS2forgbool <- F
  ###### We determine the relevant hull and the maximum 'maxpt' in this hull
  locchull <- providefullhull(notinKgspace)[[1]] ## 'full' dimensional, contains effectedConstraints
  rosglobal <- blackbox.getOption("rosglobal")
  if (length(notinKgspace)==0){ ## no extra composite hull
    maxpt <- rosglobal ## no change of hull
  } else if (length(notinKgspace)==1) {
    ## ML in composite space hull
    initpt <- fromFONKtoanyspace(rosglobal$par, colnames(locchull$vertices))
    if( ! (isPointInCHull(initpt, constraints=locchull[c("a", "b")]))) { ## rosglobal not in composite hull
      initpt <- fitobject$x[which.max(fitobject$fitted.values), ]
      initpt <- fromFONKtoanyspace(initpt, colnames(locchull$vertices))
    }
    ## regretably, optim() still stops far from the maximand even with a gradient visibly not null; fnscale and parscale have no effect
    ## testcase is g05 avec binning 2, job 72, LRT for Nb, which doesn ot find maximum when it starts from fitobject$x[which.max(fitobject$fitted.values), ] rather than from rosglobal$par
    maxpt <- optimWrapper( ##purefn, ##  parscale is provided within optimWrapper
      initval=initpt, gr=NULL,
      chullformats=locchull,
      control=list(fnscale=-1/blackbox.getOption("scalefactor"), trace=FALSE, maxit=10000)) ## returns in fittedNames space
    ## maxpt is in kg space
  } else if (length(notinKgspace)>1) {
    stop.redef("Several composite variables: case not handled in LRT computation")
  }
  ## maxpt can still be improved below so we should not yet determine Promaxval
  ###### We determine the likelihood for the tested parameters
  if(length(LRTfixedvals)==fittedparamnbr) { ## no profiling needed, generic LRT
    if (isPointInCHull(unlist(LRTfixedvals), constraints=locchull[c("a", "b")])) {
      KrigVec <- tofullKrigingspace(LRTfixedvals, fixedlist=NULL)
      profpt <- list(par=LRTfixedvals, value=purefn(KrigVec, testhull=F))
    } else {
      KrigVec <- NA
      message.redef("(!) LRT sought for a point not in the convex envelope of the sampled parameter points. No value returned.")
      profpt <- list(par=LRTfixedvals, value=NA)
    }
  } else { ## We determine the profile likelihood for the tested parameters
    ## constrained optimization in locchull determined above by names(LRTfixedvals)
    ## (this hull can be redetermined within profile(), but we can check whether the results are the same if we transmit it as argument)
    profpt <- profileBySubHull(fixedlist=LRTfixedvals, locchull=locchull, templateinkg=maxpt$par, max.only=F, usezoom=T)
    if( ! any(is.na(profpt$par)) ) {
      KrigVec <- tofullKrigingspace(profpt$par, fixedlist=LRTfixedvals)
    } else {
      KrigVec <- NA
      message.redef("(!) No profile found for given testPoint. Maybe not in convex envelope of sampled parameter points? No value returned.")
    }
    # profile() minimal return is *canonical* $par, $value, $message
    ## if higher maximum found...
    if ( ! is.na(profpt$value) && profpt$value>maxpt$value) { ## profile better than ML
      ## then we seek a higher maximum around this value (in locchull)
      ## reconversion to kriging space for input to optim(Krigvec, purefn, ...) below
      ros <- optimWrapper( ##purefn,
        initval=KrigVec, gr=NULL,
        chullformats=locchull,
        control=list( ## parscale is provided within optimWrapper
          fnscale=-1/blackbox.getOption("scalefactor"), trace=FALSE, maxit=10000))
      if (maxpt$value>ros$value) { ## we already know that profpt > maxpt and we should have obtained ros > profpt
        cat("(!) Problem: profpt > maxpt > ros !", "\n")
        cat("Are tested values within sampled range?", "\n")
      } else { ## normal case, ros > profpt > maxpt, ros becomes the new maxpt
        maxpt <- ros ## ros is automatically in the relevant hull 'locchull'
        ## and if this is relevant for rosglobal in Kriging hull...
        if ( ros$value>rosglobal$value ## if maximum lik is higher than that of current rosglobal
             &&  isPointInCHull(ros$par, constraints=blackbox.getOption("hulls")$Kgtotal[c("a", "b")]) ) { ## and maximum point is also in Kriging hull
          canonized <- canonize(ros$par)
          DemographicModel <- blackbox.getOption("DemographicModel")
          if ("IBD" %in% DemographicModel) {
            ros <- c(ros, list(latt2Ns2=canonized$latt2Ns2))
          } else {
            if("OnePopVarSize" %in% DemographicModel) {
              ros <- c(ros, list(Nratio=canonized$Nratio))
            } else 
                if ("OnePopFounderFlush" %in% DemographicModel) {ros <- c(ros, list(Nratio=canonized$Nratio), list(NactNfounderratio=canonized$NactNfounderratio), list(NfounderNancratio=canonized$NfounderNancratio))}
          }
          ros <- c(ros, canonVP=list(canonized$canonVP)) ## !!comme dans findGlobalMLE !!
          blackbox.options(rosglobal=ros)
        }
      }
    } ## endif profpt$value>maxpt$value
  }
  Promaxval <- maxpt$value
  Protest <- profpt$value
  if (length(usernames)>1) {
    namestring <- paste("(", paste(usernames, collapse=", "), ")", sep="")
    valstring <- paste("(", paste(uservalues, collapse=", "), ")", sep="")
    fixedstring <- paste(namestring, " = ", valstring, sep="")
  } else fixedstring <- paste(usernames, " = ", uservalues, sep="")
  if (is.na(Protest)) {
    message.redef(paste("From LRT computation: no profile likelihood found for value [ ", fixedstring, " ] (out of envelope of kriged points ?)"))
    LRT <- NA
  } else LRT <- -2*(Protest-Promaxval)
  if (!is.na(LRT) && LRT>0) {
    pval <- (1-pchisq(LRT, df=length(LRTfixedvals))) ## p values
  } else {pval <- NA}
  if( ! any(is.na(KrigVec))) {
    profpt.canon <- canonize(KrigVec)$canonVP
    resust <- paste("LRT for [ ", fixedstring, " ]: ", prettynum(LRT), ", and Pvalue: ", prettynum(pval), #
                  " (at [ ", paste(prettynum(profpt.canon), collapse=", "), " ]_c )", sep="")
    cat(resust, "\n")
    write(resust, file=cleanResu)
  }
  ## For the one case where rosglobal can have been improved
  returncode <- profpt$edgelevel
  if(is.null(returncode)) {
    returncode <- NA ## sinon une colonne manquera dans la sortie fichier ##
  } else if (returncode>0) {
    write("(!) Likelihood maximum for the tested value is at the edge of the parameter space used for testing this variable.", file=cleanResu)
  }
  tmp <- maxpt$edgelevel # not rosglobal$edgelevel  if maximization in non-kriging space
  if (tmp>0) {
    write("(!) Global likelihood maximum is at the edge of the parameter space used for testing this variable.", file=cleanResu)
  }
  if (tmp>0) returncode <- returncode+tmp/(10^ceiling(log(tmp, 10))) ## second summand goes in decimal part of returcode
  GOP <- blackbox.getOption("ptNbrforCI")
  writeoutput(locst, returncode=returncode, LRT, pval, GOP)
  if (condS2forgbool) {
    ## nothing to do to profpt$par (this does not contain the value at which the profile is sought!)
    ## nothing to do to maxpt$par this is always in kriging space...
  }
  return(list(LRTnames=names(LRTfixedvals), ## standardNames
              LRTfixedvals=LRTfixedvals, ## logscale
              returncode=returncode, LRT=LRT, pval=pval,
              profpt=list(par=profpt$par, value=profpt$value), maxpt=list(par=maxpt$par, value=maxpt$value), ptNbrforCI=GOP))
} ## end LRTfn
