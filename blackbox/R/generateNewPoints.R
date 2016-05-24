set_extrapol_dlr_from_LRTs <- function(locLRTlist,dlr,FONKgNames,lowextrapol=1.1,verbose=FALSE) {
  INFO <- list(hullExpandFactor=blackbox.getOption("hullExpandFactor"))
  if (length(locLRTlist)>0 ) {
    if (verbose) { message.redef("Using locLRTlist...")  }
    ## One problem here is that LRT's of Ds2 in (twoNm, condS2) space that can return (Ds2, condS2) values implying very high twoNm values
    ## We try whether we can simply exclude generically LRT's with suspiciously low p-values on variables out of kriging space
    ## FR->FR also NA pval if profpt$value>maxpt$value LR very slightly negative (actually null). Usually caught before in kriging space but not in extra space..
    notsuspect <- unlist(lapply(locLRTlist,
                                function(l) {return( ! ( is.na(l$pval) || (l$pval<1e-6 && (length(l$LRTnames %w/o% FONKgNames)>0))))}))
    locLRTlist <- locLRTlist[notsuspect]
    LRTdlr <- unlist(lapply(locLRTlist, function(l) {return(l$LRT)})) ## should not contain NA, see computation of ...LRTlist
    LRTdlr <- max(LRTdlr)/2 ## ## >0: 10/2015
    if (LRTdlr>dlr) { ## >, not < : 10/2015
      message.redef("     LRT tested value seems outside of putative training points;")
      message.redef("     training points are extended to cover LRT tested value.")
      dlr <- LRTdlr ## extend so as to cover LRT point.
      extrapol <- lowextrapol ## but goes less far beyond the already considerably extended range ##FR->FR quite ad hoc value
    } else extrapol <- INFO$hullExpandFactor ## if LRTdlr> ...
  } else extrapol <- INFO$hullExpandFactor ## no locLRTlist, but expand !
  return(list(extrapol=extrapol,dlr=dlr))
}

provideHullwFallBack <- function(fitobject,
                                 lrthreshold, ## always used; may be modified if too few high points
                                 rosglobal, ## used only if too few high points
                                 dlrtolerance, ## not modified; used only if too few high points
                                 probErr, ## not modified; used only if too few high points
                                 verbose=FALSE
) {
  fitX <- fitobject$x
  tmpy <- fitobject$fitted.values
  dessus <- which(tmpy>lrthreshold) ## defines points with sufficiently high likelihood
  ## Apply a fall-back mechanism when there are not enough such points:
  if (length(dessus)<(2*ncol(fitX))) {
    warn <- TRUE
    dessus <- order(tmpy, decreasing=TRUE)[1:(2*ncol(fitX))]
    knots <- fitX[dessus, , drop=FALSE]
    posdlr1 <- (rosglobal$value-tmpy[dessus][2*ncol(fitX)])*dlrtolerance
    lrthreshold <- ( rosglobal$value -posdlr1 -probErr )
    if(verbose) {
      message.redef("NEW dLnL threshold value (Chi-square threshold + probable prediction error): ")
      locstring <- paste(prettynum(posdlr1 +probErr), "( = ", prettynum(posdlr1),
                         " + ", prettynum(probErr), " )", sep="")
      message.redef(locstring)
    }
  } else {
    knots <- fitX[dessus, , drop=FALSE]
    warn <- FALSE
  }
  return(list(knots=knots, lrthreshold=lrthreshold,warn=warn))
}

calcKnotsInfoWrapper <- function(fitobject,rosglobal,
                                dlr,
                                dlrtolerance=1.2, ## expand beyond *predicted* dlr threshold  ## this should be disconnected from migraine's hullExpandFactor
                                lrthreshold=NULL,
                                probErr= blackbox.getOption("RMSpred")*qnorm(0.99, 0, 1), ## qnorm thus misses 1% unidirectionally
                                knots_canon=NULL,
                                locLRTlist,FONKgNames,verbose=FALSE) {
  INFO <- blackbox.options()[c("fittedNames","samplingSpace","FONKgNames","FONKgScale","samplingScale")]
  if (is.null(knots_canon)) {
    ## ****************creating 'knots'*****************
    ## the idea in the next lines of code is that the likelihood at the true maximum may be underestimated
    ## We then add a rough estimate of the 'probable error' of this single point prediction to dlr
    ## so that the current prediction of the likelihood in that point is within the range of points
    ## that would be contained within the dlrwithError range.      [qnorm(0.99...)=2.326... ]
    dlrwithError <- dlr*dlrtolerance + probErr
    if (is.null (lrthreshold)) lrthreshold <- rosglobal$value - dlrwithError
    if (verbose) {
      message.redef("Creating knots from fitobject...")
      message.redef("Default dLnL threshold value (Chi-square threshold*tolerance + probable prediction error): ")
      locstring <- paste(prettynum(dlrwithError), "( = ", prettynum(dlr),"*",prettynum(dlrtolerance), " + ", prettynum(probErr), " )", sep="")
      message.redef(locstring)
    }
    locblob <- provideHullwFallBack(fitobject=fitobject, ## selection of points in fitobject
                                    lrthreshold=lrthreshold, ## simply returns points with higher logL, if enough of them 
                                    ## fallback args:
                                    rosglobal=rosglobal, ## used to compute LR threshold; $par NOT added
                                    dlrtolerance=dlrtolerance,probErr=probErr,verbose=verbose)
    if (locblob$warn) {
      message.redef("(!)  Few good training points for generating next points.")
      message.redef("     It is advised to compute more points.")
    }
    knots <- locblob$knots
    lrthreshold <- locblob$lrthreshold ## provideHullwFallBack input value unless fall back mechanisms is used
    knots <- rbind(knots, rosglobal$par)
    refpos <- nrow(knots) ## so that order of points matters for call of symmetricKnots()
    ## following points not used by symmetricKnots():
    keypts <- knotsFromLRTlist(locLRTlist=locLRTlist, lrthreshold=lrthreshold)
    if( ! is.null(nrow(keypts))) {knots <- rbind(knots, keypts)}
    profedges <- knotsFromProfileEdges(lrthreshold=lrthreshold)
    knots <- rbind(knots, profedges)
  } else {
    if (verbose) { message.redef(paste("Creating knots from argument 'knots_canon' of dimensions",paste(dim(knots_canon),collapse="x"),"..."))  }
    knots <- fromCanonToFONK(knots_canon)
    knots <- knots[,INFO$fittedNames]
    # lrthreshold unchanged
  }
  ## knots must now be in sampling space (and scale) for comparison to constraints in LowUpFromKnots -> LowUpfn
  samplingSpace <- INFO$samplingSpace
  NOTINKGSPACESCALE <- (length(samplingSpace %w/o% INFO$FONKgNames)>0
                        || ! identical(INFO$FONKgScale,
                                       INFO$samplingScale))
  if (NOTINKGSPACESCALE) {
    knots <- t(apply(knots, 1, function(v) {
      fromFONKtoanyspace(v, samplingSpace,outputScale=INFO$samplingScale) ## FR->FR pas vrai si fixedparam => calcKntsInfo => redundant => convhulln plante sur colonne constante
    }))
  }
  LowUp <- LowUpFromKnots(knots=knots,verbose=FALSE) ## not verbose in all cases
  localLow <- LowUp$localLow
  localUp <- LowUp$localUp
  if (rosglobal$edgelevel>1) {
    if (verbose) { message.redef("Using symmetricKnots...")  }
    # takes symmetric points in sampling space as the first four args are in this space
    #### reflexion if nvecessary, of the knotsFromFitobject, not of the keypts; note that reflected points are possibly redundant
    fitobjectMirror <- constrSymOrReflWrapper(points=knots[1:(refpos-1),,drop=FALSE],
                                              logLcheck=lrthreshold, refpoint=unlist(knots[refpos,]),
                                              lower=localLow, upper=localUp) ## was reflexion of knots through mirrorface, now symmetry through refpoint...
    knots <- rbind(knots, fitobjectMirror)
  }
  knots <- shrink_knots(knots,lower=localLow, upper=localUp, verbose=verbose)
  knotsInfo <- try(calcKnotsInfo(knots, previousknotsInfo=NULL),silent=TRUE)
  if (inherits(knotsInfo,"try-error")) { 
    amknots <- as.matrix(knots) ## the sweep syntax does not work on data.frame's 
    colvecs <- t(sweep(amknots[-1,,drop=FALSE],2,amknots[1,],`-`))
    qrvecs <- qr(colvecs)
    if ((rk <- qrvecs$rank)<ncol(knots)) { ## handles lower-dimensional knots (cf explanation below end of function)
      shrink_x <- shrink_knots(fitobject$x,lower=localLow, upper=localUp, verbose=verbose) 
      colx <- t(sweep(shrink_x,2,knots[1,],`-`)) # col vectors of shrink_x relative to first knot 
      projx <- qr.qty(qrvecs,colx) ## col vectors
      ## now we need to find good points that will increase the dimension of completed knots 
      order_fity <- order(fitobject$fitted.values,decreasing = TRUE) 
      missingdims <- seq(ncol(knots)-rk) 
      ## note shrink_x ordered according to loglik unshrinked values
      blob <- cbind(shrink_x, t(projx[rk+missingdims,,drop=FALSE]))[order_fity,]
      additionalpts <- integer(0)
      for (colit in (ncol(fitobject$x)+missingdims)) {
        outlying <- which(abs(blob[,colit])>1e-6) ## setdiff to make sure a distinct point is selected for each dimension: 
        additionalpts <- c(additionalpts,setdiff(outlying,additionalpts)[1L])
      }
      knots <- rbind(knots,shrink_x[order_fity[additionalpts],])
      knots <- shrink_knots(knots,lower=localLow, upper=localUp, verbose=verbose)
      knotsInfo <- calcKnotsInfo(knots, previousknotsInfo=NULL)
    } else return(knotsInfo) ## returns the error: unhandled case
  }
  knotsInfo$lrthreshold <- lrthreshold
  knotsInfo$LowUp <- LowUp
  return(knotsInfo)
} ## end def calcKnotsInfoWrapper

if (FALSE) { ## explanation with messy code
  rowpts <- cbind(rep(1,5),c(0,1,-1,1,-1),c(0,1,1,-1,-1)) ## : 2D square in 3D coords
  ## define first point as 'origin'
  colvecs <- t(sweep(rowpts[-1,,drop=FALSE],2,rowpts[1,],`-`)) # t() -> col vectors
  qrvecs <- qr(colvecs)
  good <- seq(qrvecs$rank)
  lowerdimensionalprojections <- (t(qr.Q(qrvecs)) %*% cbind(0,colvecs))[good,] ## putting back the origin
  # in example, clearly origin + (+/-)sqrt(2) projections on obvious orthogonal basis
  indices <- convhulln(t(lowerdimensionalprojections),"Pp") ## example: eliminates the central point which is the origin
  # put back the full dimension, which first point still 'origin':
  fulldimensionalfromorigin <-  qr.Q(qrvecs)[,good] %*% lowerdimensionalprojections
  # put back the coordinates of the origin:
  colvecsagain <- (sweep(t(fulldimensionalfromorigin)[,,drop=FALSE],2,rowpts[1,],`+`))
  ## see the selected points:
  colvecsagain[unique(as.vector(resu)),] ## = points in the smallest convex simplicial complex
}

generateNewPoints <- function(pointnbr=-1, ##   pointnbr is the target size for cumulgoodpoint WITHIN this call
                              step, previous=NULL,
                              extrapol=if(step=="expand") {blackbox.getOption("hullExpandFactor")} else {1},
                              knotsInfo, ## something generated by calcKnotsInfo,with lrthreshold added, (calcKnotsInfoWrapper)
                              focalPoint=NULL,
                              fitforEI=NULL,
                              vT=volTriangulation(knotsInfo$knotsVH$vertices), 
                              verbose=FALSE ## verbosity for development purposes
                              ) {
#browser()
  INFO <- blackbox.options()[c("FONKgLow","fittedNames","fittedparamnbr","ParameterNames","FONKgNames","memcheck","barycenterFn","scalefactor",
                               "parDigits")]
  FONKgLow <- INFO$FONKgLow
  fittedNames <- INFO$fittedNames
  fittedparamnbr <- INFO$fittedparamnbr
  if (pointnbr<fittedparamnbr+1) return()
  othernames <- INFO$FONKgNames %w/o% fittedNames
  knots <- knotsInfo$knotsVH$vertices
  lrthreshold <- knotsInfo$lrthreshold
  ## *************** knots may be contracted later ************* ## but vT may not.
  blob <- LowUpFromKnots(knots=knots,verbose=verbose)
  localLow <- blob$localLow;localUp <- blob$localUp #
  canonNames <- INFO$ParameterNames
  if (is.null(previous)) {
    previous <- as.data.frame(matrix(nrow=0, ncol=length(canonNames)))
    colnames(previous) <- canonNames
  }
  nini <- NROW(previous)
  nbrtriedpoints <- 0;successrate <- 0
  nbrfeasibletriedpoints <- 0;feasiblesuccessrate <- 0
  prevmsglength <- 0
  ##To allow expansion beyond the knots, but still constrained by localLow, localUp:
  ## extrapolation occurs where knots are more constraining than localLow, localUp
  newLowUp <- FALSE   ## indicates that we will fit the knots in Low and Up
  minFeasibleSuccessrate <- 0.05
  maximumunfeasiblerate <- 0.5
  suspiciouslyLowRate <- (minFeasibleSuccessrate*maximumunfeasiblerate/2.5)
  while (nrow(previous)<pointnbr) {
    if (interactive() && RAMavail() < INFO$memcheck) { ## testing memory available for malloc
      message("RAM available < blackbox.getOption(\"memcheck\"), use the browser for diagnosis:")
    }
    trynbr <- pointnbr-nrow(previous)
    if(successrate>0) {
      trynbr <- trynbr/successrate;
      trynbr <- ceiling(trynbr+2*sqrt(trynbr)) ## floor can yield 0 and an infinite loop
      trynbr <- min(trynbr, max(2000L,pointnbr))
    }
    if(newLowUp) {  ## can occur only in "expand" step
      #print(c("avant calcKnotsInfo",RAMavail(),blackbox.getOption("memcheck")))
      knots <- shrink_knots(knots,lower=localLow, upper=localUp, verbose=verbose)
      knotsInfo <- calcKnotsInfo(knots, previousknotsInfo=knotsInfo)
      knots <- knotsInfo$knotsVH$vertices
      newLowUp <- FALSE
    } ## else leave knots unchanged
    if (is.null(focalPoint)) { ## then perform an optim to find default focalPoint
      barycenter <- do.call(INFO$barycenterFn, list(vertices=knots))
      barycenter <- tofullKrigingspace(barycenter) ## as the knots (presumably)
      optr <- optimWrapper( ##purefn,
        initval=barycenter, gr=NULL,
        chullformats=knotsInfo$knotsVH,
        control=list( ## parscale is provided within optimWrapper
          fnscale=-1/INFO$scalefactor, trace=FALSE, maxit=10000))
      ## may fail if optimize on a facet of the hull
      focalPoint <- optr$par
    }
    #}
    #
    if (step=="expand") {
      candidates <- rExpandedHull(trynbr,
                                knotshull=knots,
                                hrepr=knotsInfo$knotsVH$Hrep,
                                vT=vT,
                                extrapol=extrapol,
                                focalPoint=focalPoint)
      nbrtriedpoints <- nbrtriedpoints+nrow(candidates)
      candidatesBlob <- select_ByHull_ByLogic(candidates,
                                              localLow=localLow,
                                              localUp=localUp,
                                              fittedparamnbr=fittedparamnbr,
                                              fittedNames=fittedNames,
                                              maximumunfeasiblerate=maximumunfeasiblerate)
      candidates <- candidatesBlob$candidates
      nbrfeasibletriedpoints <- nbrfeasibletriedpoints+nrow(candidates)
      newLowUp <- candidatesBlob$newLowUp
      if (newLowUp) {
        localLow <- candidatesBlob$localLow
        localUp <- candidatesBlob$localUp
      }
      # !! selectByLR calls toCanonical => rajoute les fixed Pars
      #print(c("avant selectByLR",RAMavail(),blackbox.getOption("memcheck")))
      candidates <- selectByLR(candidates=candidates,
                               lrthreshold=lrthreshold,
                               othernames=othernames,
                               FONKgLow=FONKgLow,
                               fittedparamnbr=fittedparamnbr)
      # !! selectByLR calls toCanonical => rajoute les fixed Pars
    } else {
      if ( step=="fillbyEI" ) {
        #print(c("avant rhullByEI",trynbr, RAMavail(),blackbox.getOption("memcheck")))
        candidates <- rhullByEI(n= trynbr,
                                vT=vT, object=fitforEI)
      } else if ( step=="fill" ) {
        candidates <- rhullByvT(trynbr, vT)
      }
      nbrtriedpoints <- nbrtriedpoints+nrow(candidates)
      ## this samples within vT => does not extrapolate beyond vT
      candidates <- select_ByLogic(candidates,fittedNames=fittedNames,fittedparamnbr=fittedparamnbr)
      nbrfeasibletriedpoints <- nbrfeasibletriedpoints+nrow(candidates)
      candidates <- toCanonical(candidates, FONKgLow=FONKgLow, othernames=othernames)
    }
    #
    if (! is.null(candidates)) {
      candidates <- signif(candidates, INFO$parDigits)
      previous <- unique(rbind(previous, candidates)) ## previous accumulates valid results
      msg <- paste("    ...already", min(pointnbr, nrow(previous)), "points generated...")
      prevmsglength <- overcat(msg, prevmsglength)
    }
    n_added <- NROW(previous)-nini ## added to 'previous' since input 'previous'
    successrate <- n_added/nbrtriedpoints
    feasiblesuccessrate <- n_added/nbrfeasibletriedpoints
    if (n_added==0L || successrate<suspiciouslyLowRate) {
      if(nbrtriedpoints>(pointnbr/suspiciouslyLowRate)) { ## if successrate remains too low after the second round
        if (step!="expand") {
          message.redef("From generateNewPoints(...): too few points satisfy constraints.")
          message.redef("This may be due to a non-default krigScale setting. In that case,")
          message.redef("  use only samplingScale, which also control krigScale by default.")
          stop.redef("I exit. See messages from generateNewPoints(...).")
        } else {
          break; ## this is a reasonable outcome for expand step. We keed available points and will call another innerfill step.
        }
      }
    }
    if (step=="expand") {
      if (tofKpredict.nohull(focalPoint, fixedlist=NULL)>lrthreshold && feasiblesuccessrate<minFeasibleSuccessrate) {
        message.redef("");message.redef(paste("(info)  successrate: ", successrate, sep=""))
        message.redef("");message.redef("(info) Low success rate. Trying smaller expansion factor...")
        ## aiming to set feasiblesuccessrate at minFeasibleSuccessrate... but may provide too small extrapol as feasiblesuccessrate may depend on larger (ante)penultimate values
        locsuccessrate <- max(1/nbrfeasibletriedpoints,feasiblesuccessrate) ## avoid 0...
        extrapol <- extrapol*((locsuccessrate/minFeasibleSuccessrate)**(1/fittedparamnbr))
        message.redef(paste("(info)  New extrapolation factor: ", extrapol, sep=""))
        ## plus: we could also consider expanding extrapol when the success rate is too large
        #### plus potential problem here:
      }
      ##FR->FR alternatively we might expand if successrate is large; but expand through localLow, localUp is more difficult because fitinbounds can only shrink
    }
  } ##while
  if (nrow(previous)>pointnbr) previous <- previous[1:pointnbr, , drop=FALSE]
  return(previous)
} # end generateNewPoints

## complete points if not full dimentional before calling canonize
## input has dim. Since 2016/01/05, output always has dim
toCanonical <- function(candidates, FONKgLow, othernames) {
  INFO <- list(ParameterNames=blackbox.getOption("ParameterNames"))
  if (length(othernames)>0) {
    if (is.data.frame(candidates)) {
      candidates <- cbind(candidates,t(FONKgLow[othernames])) ## syntax OK for df et matrix et no rownames warning.
    } else {
      candidates <- cbind(candidates,FONKgLow[othernames])
      colnames(candidates)[ncol(candidates)+1L-seq_len(length(othernames))] <- othernames
    }
  }
  ## now this has the dimension of FONKgNames
  if(nrow(candidates)>1) {
    candidates <- apply(candidates, 1, function(v) {canonize(v)$canonVP}) ## transposed // expected; except if fittedparamnbr==1...
    if (length(INFO$ParameterNames)>1L) {
      candidates <- t(candidates)
    } else candidates <- matrix(candidates, ncol=1)
    ## apply loses $canonVP names and instead copies the candidates' names...
    colnames(candidates) <- INFO$ParameterNames
  } else {
    candidates <- t(canonize(candidates)$canonVP) ## t() converts to matrix with the column names
  }
  return(candidates)
}


selectByLR <- function(candidates, lrthreshold, othernames, FONKgLow, fittedparamnbr) {
  pred <- apply(candidates, 1, purefn, testhull=F) ## must be ins .blackbox.data$options$KrigSpace and Scale
  bons <- ( (!is.na(pred)) & (pred> lrthreshold )) ## selection according to predicted likelihood
  if(any(bons)) {
    candidates <- candidates[which(bons), , drop=FALSE]
    candidates <- toCanonical(candidates=candidates, FONKgLow=FONKgLow, othernames=othernames)
    return(candidates)
  } else return(NULL)
}

select_ByLogic <- function(candidates, fittedNames,fittedparamnbr) {
  #   candidatesinhull <- apply(candidates, 1, fitinbounds, lower=localLow, upper=localUp) ## in bounds provided by LowUpfn
  #   if (fittedparamnbr>1) {
  #     candidatesinhull <- t(candidatesinhull) ## occurs in the one-dimensional case
  #   } else candidatesinhull <- matrix(candidatesinhull, ncol=1)
  #   colnames(candidatesinhull) <- colnames(candidates) ## may be only a subset of FONKgNames, or out of fullKrigspace
  candidates <- apply(candidates, 1, tofullKrigingspace) ## back to Kriging space (but this generate NaN's)
  if (fittedparamnbr>1) {
    candidates <- t(candidates)
  } else candidates <- matrix(candidates, ncol=1)
  colnames(candidates) <- fittedNames ## now in full kriging space but not necessarily full FONKgNames
  unfeasible <- apply(candidates, 1, function(v) {any(is.nan(v))})
  candidates <- candidates[!unfeasible, , drop=FALSE]
  return(candidates)
}


select_ByHull_ByLogic <- function(candidates, localLow, localUp, fittedparamnbr, fittedNames, maximumunfeasiblerate) {
  candidatesinhull <- apply(candidates, 1, fitinbounds, lower=localLow, upper=localUp) ## in bounds provided by LowUpfn
  if (fittedparamnbr>1) {
    candidatesinhull <- t(candidatesinhull) ## occurs in the one-dimensional case
  } else candidatesinhull <- matrix(candidatesinhull, ncol=1)
  colnames(candidatesinhull) <- colnames(candidates) ## may be only a subset of FONKgNames, or out of fullKrigspace
  candidates <- apply(candidatesinhull, 1, tofullKrigingspace) ## back to Kriging space (but this generate NaN's)
  if (fittedparamnbr>1) {
    candidates <- t(candidates) ## occurs in the one-dimensional case
  } else candidates <- matrix(candidates, ncol=1)
  colnames(candidates) <- fittedNames ## now in full kriging space but not necessarily full FONKgNames
  unfeasible <- apply(candidates, 1, function(v) {any(is.nan(v))})
  newLowUp <- FALSE
  if((length(which(unfeasible))/nrow(candidates))>maximumunfeasiblerate) { ## heavy burden of unfeasible values
    ## identify the problem in hull space
    for (st in colnames(candidatesinhull)) {
      numnans <- candidatesinhull[, st]
      maxnum <- max(numnans[!unfeasible]) ## the max st value for which some meaningful points was generated
      alwaysbad <- numnans[which(numnans>maxnum)] ## extreme st values beyond the above value: a subset of NaN generating values
      if (length(alwaysbad)>0) {
        localUp[st] <- min(alwaysbad) ## will affect recomputation of knots is successrate too low
        newLowUp <- TRUE
        message.redef(paste("(info) New localUp[", st, "] value: ", localUp[st], sep=""))
      } ## else nothing can be done on the block constraints
      minnum <- min(numnans[!unfeasible]) ## the max st value for which some meaningful points was generated
      alwaysbad <- numnans[which(numnans<minnum)] ## extreme st values beyond the above value: a subset of NaN generating values
      if (length(alwaysbad)>0) {
        localLow[st] <- max(alwaysbad) ## will affect recomputation of knots is successrate too low
        newLowUp <- TRUE
        message.redef(paste("(info) New localLow[", st, "] value: ", localLow[st], sep=""))
      } ## else nothing can be done on the block constraints
    }
  }
  ## select candidates
  candidates <- candidates[!unfeasible, , drop=FALSE] ## moved outside the test, 10/2015
  return(list(candidates=candidates,
              newLowUp=newLowUp, localLow=localLow, localUp=localUp))
}

shrink_knots <- function(knots, lower, upper,verbose=FALSE) {
  if (verbose) { message.redef(paste("Shrinking the (",paste(dim(knots),collapse="x"),") knots to the bounds..."))  }
  redundantknots <- apply(knots, 1, fitinbounds, lower=lower, upper=upper) ## tentatively SHRINKING the HULL to new bounds
  if ( ! is.null(redundantknots)) {
    if (is.matrix(redundantknots)) {
      redundantknots <- t(redundantknots) ## classic apply() problem
    } else dim(redundantknots) <- c(length(redundantknots),1L) # equiv to conversion to column...
    colnames(redundantknots) <- colnames(knots)
  }
  return(redundantknots)
}

calcKnotsInfo <- function(redundantknots, previousknotsInfo) {
  INFO <- list(`redundant.mode`=blackbox.getOption("redundant.mode"))
  if ( ! is.null(previousknotsInfo) && identical(redundantknots, previousknotsInfo$redundantknots) ) {
    return(previousknotsInfo) ## unchanged
  } else {## THEN calling redundant
    redmode <- switch(INFO$`redundant.mode`,
                      "noElim"="no.elim",
                      "alwaysRational"="rational",
                      "alwaysDouble"="double",
                      "defaultPrecision"="rational")
    knotsVH <- resetCHull(redundantknots, formats=c("vertices", "constraints"), redundant.mode=redmode)
    if (inherits(knotsVH,"try-error")) {
      return(knotsVH)
    } else return(list(knotsVH=knotsVH, redundantknots=redundantknots))
  }
}
