
generateNextpointsfromCI <- function(n, ## number of points generated PER CI point
                                     CIpointsList = blackbox.getOption("CIpointsList"),
                                     posdlr, ## >0 10/2015
                                     previous=NULL, ## in canonical space
                                     fitforEI,
                                     dlrtolerance=1.2, ## expand beyond *predicted* dlr threshold ## not controlled by Migraine
                                     LowUp,
                                     verbose=FALSE
                                     ) {
  if (verbose) { message.redef("Generating points from confidence interval bounds...")  }
  fitobject <- blackbox.getOption("fitobject")
  rosglobal <- blackbox.getOption("rosglobal")
  fittedNames <- blackbox.getOption("fittedNames")
  fittedparamnbr <- blackbox.getOption("fittedparamnbr")
  FONKgNames <- blackbox.getOption("FONKgNames")
  FONKgLow <- blackbox.getOption("FONKgLow")
  othernames <- FONKgNames %w/o% fittedNames
  #
  probErr <- blackbox.getOption("RMSpred")*qnorm(0.99, 0, 1) ## qnorm thus misses 1% unidirectionally
  locblob <- provideHullwFallBack(fitobject=fitobject,
                                  lrthreshold=( rosglobal$value -posdlr*dlrtolerance -probErr ),
                                  rosglobal=rosglobal,
                                  dlrtolerance=dlrtolerance,probErr=probErr,verbose=verbose)
  if (locblob$warn) {
    message.redef("(!)  Few good training points for ascertaining CIs.")
    message.redef("     It is advised to compute more points.")
  }
  knotsFromFitobject <- locblob$knots
  #
  cumulgoodpoints <- as.data.frame(matrix(nrow=0, ncol=fittedparamnbr))
  parnames <- intersect(names(CIpointsList), fittedNames) # the following code works without this, but
  #... the current problem is that fixedPar=latt2Ns2 not in Kgspace can generate points far from the main cluster of points => ugly graphics
  for (fixedPar in parnames) {
    for (it in seq_len(nrow(stcipoints <- CIpointsList[[fixedPar]]))) {
      CIpoint <- stcipoints[it,]
      if (length(CIpoint)==1L) names(CIpoint) <- colnames(stcipoints)
      obspred <- predict(fitforEI, newdata=CIpoint, variances=list(linPred=TRUE))
      obsSE <- attr(obspred, "predVar")
      obsSE[obsSE<0] <- 0
      ## attention le signe de la prediction est inversé...
      fitforEI$Qmax <- max( - obspred+1.96 * sqrt(obsSE)) ## best improvement function for the currently predicted CI point
      ## CIpoint is is fullKrigingspace
      if ( ! fixedPar %in% names(CIpoint)) {
        fixed <- as.list(fromFONKtoanyspace(CIpoint, outputnames=fixedPar))
      } else {
        fixed <- as.list(CIpoint[fixedPar])
      }
      locvertices <-  provideVertices(fixedPar, knotsFromFitobject)
      if (ncol(locvertices)>1L) {
        subHull <- subHullWrapper(vertices=locvertices, equality=fixed) ## both arguments 'logscale' if relevant
        if (is.null(subHull$vertices)) {
          locvertices <-  providefullhull(fixedPar)[[1]]$vertices
          subHull <- subHullWrapper(vertices=locvertices, equality=fixed) ## both arguments 'logscale' if relevant
        }
        if (!is.null(subHull$vertices)) {
          ## FR->FR problem ligne suiv si les $subvertices sont distincts selon unique ou ULI mais tres proches: qhull precision error... simplex not convex.
          subvT <- volTriangulation(as.matrix(subHull$vertices))
          trypoints <- rhullByEI(n=n, # (tryn =default tryn)
                                 vT=subvT, object=fitforEI, fixed=fixed)  ## mais le Qmax utilisé est celui de l'objet global...
          ## -> utilise probade dépasser le ML actuel...
          ##: by default returns points in fittedNames space
          cumulgoodpoints <- rbind(cumulgoodpoints, trypoints)
        }
      } else cumulgoodpoints <-  rbind(cumulgoodpoints, unlist(fixed)) ## add 1DCI bound to nextpoints
    }
  }
  cumulgoodpoints <- shrink_knots(cumulgoodpoints,lower=LowUp$localLow,upper=LowUp$localUp,verbose=verbose) 
  cumulgoodpoints <- toCanonical(cumulgoodpoints, FONKgLow=FONKgLow, othernames=othernames)
  cumulgoodpoints <- unique(rbind(previous, cumulgoodpoints))
  msg <- paste("    ...already", nrow(cumulgoodpoints), "points generated...")
  cat(msg)
  return(list(cumulgoodpoints=cumulgoodpoints, knotsInfo=NULL))
}

